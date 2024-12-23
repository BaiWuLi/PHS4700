function [coup, vbf, ti, x, y, z] = Devoir2(option, rbi, vbi, wbi)
    assert(ismember(option, [1, 2, 3]), 'Option doit être 1, 2 ou 3');
    assert(norm(vbi) <= 35, 'La vitesse initiale du centre de masse de la balle ne peut jamais excéder 35 m/s');
    assert(norm(wbi) <= 940, 'La vitesse angulaire de la balle autour de son centre de masse ne peut excéder 940 rad/s');

    global Lt lt ht lf df hf mb Rb Ab p Cv CM;

    % Table
    ht = 76e-2; % hauteur de la table
    Lt = 2.74; % longueur de la table
    lt = 1.525; % largeur de la table

    % Filet
    hf = 15.25e-2; % hauteur du filet
    lf = 1.83; % largeur du filet
    df = 15.25e-2; % debordement du filet

    % Balle
    mb = 2.74e-3; % masse de la balle
    Rb = 1.99e-2; % rayon de la balle
    Ab = pi * Rb^2; % aire efficace de la balle

    % Forces
    p = 1.2; % masse volumique de l'air
    Cv = 0.5; % coefficient de frottement visqueux
    CM = 0.29; % coefficient de Magnus

    qi = q(rbi, vbi, wbi);
    ti = 0;
    dt = 0.0001;
    coup = -1;

    [coup, vbf, ti, x, y, z] = simulation(qi, ti, dt, option, coup);
    [ti, x, y, z] = compresser(ti, x, y, z);
end

function [coup, vbf, ti, x, y, z] = simulation(qi, ti, dt, option, coup)
    while coup == -1
        ti = [ti, ti(end) + dt];
        q = RungeKutta4(qi(end,:), ti(end), dt, @g, option);
        qi = [qi; q];

        coup = Coup(qi(end-1,:), qi(end,:));
    end

    if ~all(erreurs(qi) < 1e-9)
        ti(end) = [];
        qi(end,:) = [];
        dt = dt / 10;
        coup = -1;

        [coup, vbf, ti, x, y, z] = simulation(qi, ti, dt, option, coup);
        return;
    end

    qi(end,:) = [];
    vbf = qi(end, 1:3)';
    x = qi(:, 4)';
    y = qi(:, 5)';
    z = qi(:, 6)';
end

function q0 = q(rbi, vbi, wbi)
    q0 = zeros(1, 18);

    q0(1:3) = vbi;
    q0(4:6) = rbi;
    q0(7:9) = wbi;
    q0(10) = 1; % Rxx
    q0(14) = 1; % Ryy
    q0(18) = 1; % Rzz
end

function q1 = RungeKutta4(q0, t, dt, g, option)
    k1 = g(q0, t, option);
    k2 = g(q0 + dt/2*k1, t + dt/2, option);
    k3 = g(q0 + dt/2*k2, t + dt/2, option);
    k4 = g(q0 + dt*k3, t + dt, option);

    q1 = q0 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function g_ = g(q, t, option)
    global mb Rb Ab p Cv CM;

    Fg = zeros(1, 3); % Force gravitationnelle
    Fv = zeros(1, 3); % Force de frottement visqueux
    FM = zeros(1, 3); % Force de Magnus

    if option >= 1
        Fg = mb * [0, 0, -9.8];
    end
    if option >= 2
        Fv = -0.5 * p * Cv * Ab * norm(q(1:3)) * q(1:3);
    end
    if option >= 3
        FM = 4 * pi * Rb^3 * p * CM * cross(q(7:9), q(1:3));
    end

    F = Fg + Fv + FM;

    g_ = zeros(1, 18);
    g_(1:3) = F / mb; % a
    g_(4:6) = q(1:3); % v
    g_(7:9) = zeros(1, 3); % @
    g_(10) = q(8) * q(16) - q(9) * q(13); % Rxx
    g_(11) = q(8) * q(17) - q(9) * q(14); % Rxy
    g_(12) = q(8) * q(18) - q(9) * q(15); % Rxz
    g_(13) = q(9) * q(10) - q(7) * q(16); % Ryx
    g_(14) = q(9) * q(11) - q(7) * q(17); % Ryy
    g_(15) = q(9) * q(12) - q(7) * q(18); % Ryz
    g_(16) = q(7) * q(13) - q(8) * q(10); % Rzx
    g_(17) = q(7) * q(14) - q(8) * q(11); % Rzy
    g_(18) = q(7) * q(15) - q(8) * q(12); % Rzz
end

function coup = Coup(q0, q1)
    global Lt lt ht hf df Rb;

    x0 = q0(4);
    x1 = q1(4);
    y1 = q1(5);
    z0 = q0(6);
    z1 = q1(6);

    traverse_plan_filet = (x0 < (Lt/2-Rb) && (Lt/2-Rb) <= x1) || (x1 <= (Lt/2+Rb) && (Lt/2+Rb) < x0); % traverse par la gauche ou la droite
    dans_zone_filet = ((-df-Rb) < y1 && y1 < (lt+df+Rb)) && ((ht+Rb) <= z1 && z1 < (ht+hf+Rb));
    if traverse_plan_filet && dans_zone_filet
        coup = 2;
        return;
    end

    traverse_plan_table = (z0 > (ht+Rb) && (ht+Rb) >= z1);
    dans_zone_table = (0 <= x1 && x1 <= Lt) && (0 <= y1 && y1 <= lt);
    if traverse_plan_table && dans_zone_table
        dans_zone_adverse = xor(x1 > Lt/2, x1 < x0);
        if dans_zone_adverse
            coup = 0;
        else
            coup = 1;
        end
        return;
    end

    traverse_plan_sol = (z0 > Rb && Rb >= z1);
    if traverse_plan_sol
        coup = 3;
        return;
    end

    coup = -1;
end

function [dx, dy, dz] = erreurs(qi)
    global Lt ht Rb;

    q0 = qi(end-1,:);
    q1 = qi(end,:);

    coup = Coup(q0, q1);

    x = q0(4);
    z = q0(6);

    dx = abs(q1(4) - q0(4)); % erreur en x
    dy = abs(q1(5) - q0(5)); % erreur en y
    dz = abs(q1(6) - q0(6)); % erreur en z

    if coup == 0 || coup == 1
        dz = abs((ht + Rb) - z);
    elseif coup == 2 && q1(4) < Lt/2
        dx = abs((Lt/2 - Rb) - x);
    elseif coup == 2 && q1(4) >= Lt/2
        dx = abs((Lt/2 + Rb) - x);
    elseif coup == 3
        dz = abs(Rb - z);
    end
end

function [ti, x, y, z] = compresser(ti, x, y, z)
    n = 0;
    for i = 2:length(ti)
        dt = round(ti(i) - ti(i-1), 4);
        if dt ~= 0.0001
            n = i;
            break;
        end
    end

    ti = [ti(1:10:n), ti(end)];
    x = [x(1:10:n), x(end)];
    y = [y(1:10:n), y(end)];
    z = [z(1:10:n), z(end)];
end
