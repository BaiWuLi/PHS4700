function [Resultat, blocf, ballef, Post] = Devoir3(vbloci, avbloci, tl, vballei)
    global m_bloc A_bloc R_bloc_min R_bloc_max m_balle R_balle e_c;

    m_bloc = 1.2;
    A_bloc = 0.08;
    R_bloc_min = 0.5 * A_bloc * sqrt(3); % rayon de la sphere englobante
    R_bloc_max = 0.5 * A_bloc; % rayon de la sphere englobée
    r_bloc_i = [3; 3; 1];

    m_balle = 0.05;
    R_balle = 0.02;
    r_balle_i = [0; 0; 2];
    w_balle_i = [0; 0; 0];

    e_c = 0.8;

    q_bloc_i = qi(r_bloc_i, vbloci, avbloci);
    q_balle_i = qi(r_balle_i, vballei, w_balle_i);

    Q_bloc = q_bloc_i';
    Q_balle = q_balle_i';

    t = 0;
    dt = 0.0001;

    [Resultat, blocf, ballef, Post] = simulation(Q_bloc, Q_balle, t, dt, tl, -2);
end

function q = qi(ri, vi, wi)
    q = zeros(18, 1);

    q(1:3) = vi;
    q(4:6) = ri;
    q(7:9) = wi;
    q(10) = 1; % Rxx
    q(14) = 1; % Ryy
    q(18) = 1; % Rzz
end

function [Resultat, blocf, ballef, Post] = simulation(Q_bloc, Q_balle, t, dt, tl, Resultat)
    while Resultat == -2 % tant qu'il n'y a pas de collision
        t = [t; t(end) + dt];

        q_bloc = Q_bloc(end,:)';
        q_balle = Q_balle(end,:)';

        q_bloc = RungeKutta4(q_bloc, dt);
        if t(end) >= tl
            q_balle = RungeKutta4(q_balle, dt);
        end

        Q_bloc = [Q_bloc; q_bloc'];
        Q_balle = [Q_balle; q_balle'];

        Resultat = verifierCollision(q_bloc, q_balle);
    end

    if max([erreurs(Q_bloc); erreurs(Q_balle)]) > 1e-9 % si les erreurs sont trop grandes
        t(end) = [];
        Q_bloc(end,:) = [];
        Q_balle(end,:) = [];
        dt = dt / 10;

        [Resultat, blocf, ballef, Post] = simulation(Q_bloc, Q_balle, t, dt, tl, -2);
        return;
    end

    q_bloc_f = Q_bloc(end,:)';
    q_balle_f = Q_balle(end,:)';
    [v_bloc_f, w_bloc_f, v_balle_f, w_balle_f] = postCollision(q_bloc_f, q_balle_f, Resultat);

    blocf = zeros(6, 2);
    blocf(1:3, 1) = q_bloc_f(1:3);
    blocf(4:6, 1) = q_bloc_f(7:9);
    blocf(1:3, 2) = v_bloc_f;
    blocf(4:6, 2) = w_bloc_f;

    ballef = zeros(6, 2);
    ballef(1:3, 1) = q_balle_f(1:3);
    ballef(4:6, 1) = q_balle_f(7:9);
    ballef(1:3, 2) = v_balle_f;
    ballef(4:6, 2) = w_balle_f;

    Post = zeros(7, length(t));
    Post(1,:) = t';
    Post(2:4,:) = Q_bloc(:,4:6)';
    Post(5:7,:) = Q_balle(:,4:6)';
end

function q1 = RungeKutta4(q0, dt)
    k1 = g(q0);
    k2 = g(q0 + dt/2*k1);
    k3 = g(q0 + dt/2*k2);
    k4 = g(q0 + dt*k3);

    q1 = q0 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function Resultat = verifierCollision(q_bloc, q_balle)
    global R_balle R_bloc_min R_bloc_max;

    r_bloc = q_bloc(4:6);
    r_balle = q_balle(4:6);

    % collision de la balle avec le sol
    if r_balle(3) - R_balle <= 0
        Resultat = -1
        return;
    end

    % conditions de collision nécessaires
    collision_sol_possible = r_bloc(3) - R_bloc_min <= 0;
    collision_balle_possible = norm(r_bloc - r_balle) <= R_bloc_min + R_balle;

    if ~collision_sol_possible && ~collision_balle_possible
        Resultat = -2;
        return;
    end

    % conditions de collision suffisantes
    collision_sol_certaine = r_bloc(3) - R_bloc_max <= 0;
    collision_balle_certaine = norm(r_bloc - r_balle) <= R_bloc_max + R_balle;

    if collision_sol_certaine
        Resultat = 1;
        return;
    elseif collision_balle_certaine
        Resultat = 0
        return;
    end

    R = reshape(q_bloc(10:18), 3, 3);

    % collision avec le sol
    coins_bloc = coinsBloc(R, r_bloc);

    if min(coins_bloc(:, 3)) <= 0
        Resultat = 1;
        return;
    end

    % collision avec la balle
    r_balle = R' * (r_balle - r_bloc);
    point_contact = pointContact(r_balle);

    if isempty(point_contact)
        Resultat = -2;
    else
        Resultat = 0;
    end
end

function [dx, dy, dz] = erreurs(Q)
    q0 = Q(end-1,:);
    q1 = Q(end,:);

    dx = abs(q1(4) - q0(4));
    dy = abs(q1(5) - q0(5));
    dz = abs(q1(6) - q0(6));
end

function [v_bloc_f, w_bloc_f, v_balle_f, w_balle_f] = postCollision(q_bloc_f, q_balle_f, Resultat)
    global m_balle R_balle m_bloc A_bloc e_c;

    if abs(Resultat) == 1
        v_bloc_f = q_bloc_f(1:3);
        w_bloc_f = q_bloc_f(7:9);
        v_balle_f = q_balle_f(1:3);
        w_balle_f = zeros(3, 1);
        return;
    end

    R = reshape(q_bloc_f(10:18), 3, 3);

    % La balle est l'objet a
    r_balle = R' * (q_balle_f(4:6) - q_bloc_f(4:6));
    v_balle = R' * q_balle_f(1:3);

    % Le bloc est l'objet b
    v_bloc = R' * q_bloc_f(1:3);
    w_bloc = R' * q_bloc_f(7:9);

    r_contact = pointContact(r_balle);
    r_balle_contact = r_contact - r_balle;
    r_bloc_contact = r_contact;

    normale = - r_balle_contact / norm(r_balle_contact);
    v_balle_contact = v_balle;
    v_bloc_contact = v_bloc + cross(w_bloc, r_bloc_contact);
    v_rel_normale = dot(v_balle_contact - v_bloc_contact, normale);

    I_bloc = 1/6 * m_bloc * A_bloc^2 * eye(3);
    G_bloc = G(normale, I_bloc, r_bloc_contact);
    alpha = 1/m_balle + 1/m_bloc + G_bloc;
    alpha = 1/alpha;

    J = - alpha * (1 + e_c) * v_rel_normale * normale;

    v_balle_f = v_balle + J / m_balle;
    v_bloc_f = v_bloc - J / m_bloc;
    w_bloc_f = w_bloc - inv(I_bloc) * cross(r_bloc_contact, J);

    v_balle_f = R * v_balle_f;
    v_bloc_f = R * v_bloc_f;
    w_bloc_f = R * w_bloc_f;
    w_balle_f = zeros(3, 1);
end

function dqdt = g(q)
    dqdt = zeros(18, 1);

    dqdt(1:3) = [0; 0; -9.8]; % accélération
    dqdt(4:6) = q(1:3); % vitesse
    dqdt(7:9) = [0; 0; 0]; % vitesse angulaire

    R = reshape(q(10:18), 3, 3); % matrice de rotation
    W = [0, -q(9), q(8); q(9), 0, -q(7); -q(8), q(7), 0]; % matrice de vitesse angulaire
    dRdt = R * W;
    dqdt(10:18) = dRdt(:); % dérivée de la matrice de rotation
end

function coins_bloc = coinsBloc(rotation, translation)
    global A_bloc;
    coins_bloc = [];
    directions = [-1, 1] * A_bloc / 2;

    for i = directions
        for j = directions
            for k = directions
                coin_bloc = [i; j; k];
                coin_bloc = rotation * coin_bloc + translation;
                coins_bloc = [coins_bloc; coin_bloc'];
            end
        end
    end
end

function point_contact = pointContact(r_balle)
    global A_bloc R_balle;
    points_contact_possibles = [];

    for i = -1:1
        for j = -1:1
            for k = -1:1
                if i == 0 && j == 0 && k == 0
                    continue;
                end
                points_contact_possibles = [points_contact_possibles; i, j, k];
            end
        end
    end

    points_contact_possibles = points_contact_possibles * A_bloc / 2;

    for i = 1:size(points_contact_possibles, 1)
        point = points_contact_possibles(i, :)';
        distance_au_carre = 0;
        intervalle_valide = true;

        for j = 1:3
            if abs(point(j)) == 0.5 * A_bloc
                distance_au_carre = distance_au_carre + (r_balle(j) - point(j))^2;
            else
                intervalle_valide = intervalle_valide && abs(r_balle(j)) < A_bloc / 2;
                point(j) = r_balle(j);
            end
        end

        if intervalle_valide && distance_au_carre <= R_balle^2
            point_contact = point;
            return;
        end
    end

    point_contact = [];
end

function resultat = G(n, I, r)
    resultat = dot(n, cross(inv(I) * cross(r, n), r));
end
