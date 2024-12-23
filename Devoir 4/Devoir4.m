function [xi, yi, zi, face] = Devoir4(nout, nin, poso)
    global ELLIPSOIDE FACE1 FACE2 FACE3 FACE4 FACE5 FACE6 CONTRAINTES1 CONTRAINTES2 CONTRAINTES3 CONTRAINTES4 CONTRAINTES5 CONTRAINTES6 M N MAX_REFLEXIONS;

    % Ellipsoide = [x0, y0, z0, a, b, c]
    ELLIPSOIDE = [4, 4, 11, 3, 3, 9];

    % Plan = [a, b, c, d]
    FACE1 = [1, 0, 0, -3];
    FACE2 = [1, 0, 0, -4];
    FACE3 = [0, 1, 0, -3];
    FACE4 = [0, 1, 0, -5];
    FACE5 = [0, 0, 1, -12];
    FACE6 = [0, 0, 1, -17];

    % Limite du plan = [x_min, x_max; y_min, y_max; z_min, z_max]
    CONTRAINTES1 = [3, 3; 3, 5; 12, 17];
    CONTRAINTES2 = [4, 4; 3, 5; 12, 17];
    CONTRAINTES3 = [3, 4; 3, 3; 12, 17];
    CONTRAINTES4 = [3, 4; 4, 4; 12, 17];
    CONTRAINTES5 = [3, 4; 3, 5; 12, 12];
    CONTRAINTES6 = [3, 4; 3, 5; 17, 17];

    M = 100; % nb angles azimutaux
    N = 100; % nb angles polaires

    MAX_REFLEXIONS = 100;

    xi = [];
    yi = [];
    zi = [];
    face = [];

    poso = poso';
    dir_rayons = directionRayons(poso);

    for i = 1:length(dir_rayons)
        dir = dir_rayons(i, :)';
        [pos_image, face_atteinte] = tracerRayon(poso, dir, nin, nout);

        if isempty(pos_image)
            continue;
        end

        xi(end+1) = pos_image(1);
        yi(end+1) = pos_image(2);
        zi(end+1) = pos_image(3);
        face(end+1) = face_atteinte;
    end
end

function [dir_rayons] = directionRayons(poso)
    global ELLIPSOIDE M N;

    % Centre de l'ellipsoide par rapport a l'observateur
    x0 = ELLIPSOIDE(1) - poso(1);
    y0 = ELLIPSOIDE(2) - poso(2);
    z0 = ELLIPSOIDE(3) - poso(3);

    % Trouver les angles azimutaux limites en resolvant le probleme en 2D sur le plan xy
    [phi_min, phi_max] = anglesLimites(x0, y0, ELLIPSOIDE(4), ELLIPSOIDE(5));

    % Trouver les angles polaires limites en resolvant le probleme en 2D sur le plan wz
    w = [x0; y0; 0];
    w0 = norm(w);
    [theta_max, theta_min] = anglesLimites(w0, z0, ELLIPSOIDE(4), ELLIPSOIDE(6));
    theta_max = pi/2 - theta_max;
    theta_min = pi/2 - theta_min;

    dir_rayons = [];

    for i = 1:M
        phi = phi_min + (phi_max - phi_min) * (2*i - 1) / (2*M);
        for j = 1:N
            theta = theta_min + (theta_max - theta_min) * (2*j - 1) / (2*N);
            dir_rayons = [dir_rayons; [sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)]];
        end
    end
end

function [pos_image, face_atteinte] = tracerRayon(poso, dir_ray, nin, nout)
    global MAX_REFLEXIONS;

    dir_image = dir_ray;
    pos_image = [];
    face_atteinte = false;

    distance = distanceIntersectionEllipsoide(poso, dir_ray);
    if distance == Inf
        return;
    end

    point_incidence = poso + distance * dir_ray;
    intersection_bloc = false;

    for i = -1:MAX_REFLEXIONS % -1, car la 1ere reflexion interne a lieu apres le 2eme rebond
        dir_ray = directionReflexionRefraction(point_incidence, dir_ray, nin, nout);

        t = distancesIntersections(point_incidence, dir_ray);

        [t_min, t_argmin] = min(t); % t_argmin est le numero de la face intersectee

        if t_min == Inf % Pas d'intersection, le rayon sort de l'ellipsoide
            return;
        end

        point_incidence = point_incidence + t_min * dir_ray;
        intersection_bloc = t_argmin <= 6;
        distance = distance + t_min;

        if intersection_bloc
            pos_image = poso + distance * dir_image;
            face_atteinte = t_argmin;
            return;
        end
    end
end

function [angle_min, angle_max] = anglesLimites(x0, y0, a, b)
    % (x0, y0) est le centre de l'ellipse
    % a, b sont les demi-axes de l'ellipse
    % Equation de l'ellipse: (x - x0)^2/a^2 + (y - y0)^2/b^2 = 1
    % Equation de la droite tangent a l'ellipse passant par l'origine: y = m*x
    % On plug y = m*x dans l'equation de l'ellipse et on obtient une equation du 2eme degre en m
    A = x0^2 - a^2;
    B = -2*x0*y0;
    C = y0^2 - b^2;
    DELTA = B^2 - 4*A*C;

    assert(DELTA >= 0, "L'observateur est a l'interieur de l'ellipse");

    m1 = (-B - sqrt(DELTA)) / (2*A); % pente de la tangente inferieure
    m2 = (-B + sqrt(DELTA)) / (2*A); % pente de la tangente superieure

    angle_min = atan(m1);
    angle_max = atan(m2);
end

function t_ellipsoide = distanceIntersectionEllipsoide(pos, dir)
    global ELLIPSOIDE;

    [x0, y0, z0] = deal(pos(1), pos(2), pos(3));
    [vx, vy, vz] = deal(dir(1), dir(2), dir(3));
    [x1, y1, z1] = deal(ELLIPSOIDE(1), ELLIPSOIDE(2), ELLIPSOIDE(3));
    [a, b, c] = deal(ELLIPSOIDE(4), ELLIPSOIDE(5), ELLIPSOIDE(6));

    % Pour trouver l'intersection entre le rayon et l'ellipsoide,
    % on plug l'equation parametrique du rayon dans l'equation de l'ellipsoide :
    % ((x0 + t*vx) - x1)^2/a^2 + ((y0 + t*vy) - y1)^2/b^2 + ((z0 + t*vz) - z1)^2/c^2 = 1
    % On obtient une equation du 2eme degre en t où
    A = vx^2/a^2 + vy^2/b^2 + vz^2/c^2;
    B = 2 * ((x0 - x1)*vx/a^2 + (y0 - y1)*vy/b^2 + (z0 - z1)*vz/c^2);
    C = ((x0 - x1)/a)^2 + ((y0 - y1)/b)^2 + ((z0 - z1)/c)^2 - 1;
    DELTA = B^2 - 4*A*C;

    if DELTA <= 0
        t_ellipsoide = Inf;
        return;
    end

    t1 = (-B - sqrt(DELTA)) / (2*A);
    t2 = (-B + sqrt(DELTA)) / (2*A);
    t = sort([t1, t2]);
    t = t(t > 1e-9); % On ne garde que les t > 0

    if isempty(t)
        t_ellipsoide = Inf;
    else
        t_ellipsoide = min(t);
    end
end

function dir_rayon = directionReflexionRefraction(pos, dir, nin, nout)
    normal_i = normaleEllipsoide(pos);

    if dot(dir, normal_i) > 0 % Le rayon est a l'interieur de l'ellipsoide
        normal_i = -normal_i;
    else % Le rayon est a l'exterieur de l'ellipsoide
        [nin, nout] = deal(nout, nin);
    end

    tangente_j = cross(dir, normal_i);
    tangente_j = tangente_j / norm(tangente_j);
    tangente_k = cross(normal_i, tangente_j);

    sin_angle_incidence = dot(dir, tangente_k);
    sin_angle_refraction = (nin / nout) * sin_angle_incidence;

    if sin_angle_refraction > 1 % Reflexion interne totale
        dir_rayon = dir - 2 * dot(dir, normal_i) * normal_i;
        return;
    end

    R = coefficientReflexion(asin(sin_angle_incidence), asin(sin_angle_refraction), nin, nout);
    reflexion = rand() < R;

    if reflexion
        dir_rayon = dir - 2 * dot(dir, normal_i) * normal_i;
    else
        dir_rayon = -sqrt(1 - sin_angle_refraction^2) * normal_i + sin_angle_refraction * tangente_k;
    end
end

function t = distancesIntersections(pos, dir)
    global FACE1 FACE2 FACE3 FACE4 FACE5 FACE6 CONTRAINTES1 CONTRAINTES2 CONTRAINTES3 CONTRAINTES4 CONTRAINTES5 CONTRAINTES6;

    t_face1 = distanceIntersectionFace(pos, dir, FACE1, CONTRAINTES1);
    t_face2 = distanceIntersectionFace(pos, dir, FACE2, CONTRAINTES2);
    t_face3 = distanceIntersectionFace(pos, dir, FACE3, CONTRAINTES3);
    t_face4 = distanceIntersectionFace(pos, dir, FACE4, CONTRAINTES4);
    t_face5 = distanceIntersectionFace(pos, dir, FACE5, CONTRAINTES5);
    t_face6 = distanceIntersectionFace(pos, dir, FACE6, CONTRAINTES6);
    t_ellipsoide = distanceIntersectionEllipsoide(pos, dir);

    t = [t_face1, t_face2, t_face3, t_face4, t_face5, t_face6, t_ellipsoide];
end

function normale = normaleEllipsoide(pos)
    global ELLIPSOIDE;

    [x, y, z] = deal(pos(1), pos(2), pos(3));
    [x0, y0, z0] = deal(ELLIPSOIDE(1), ELLIPSOIDE(2), ELLIPSOIDE(3));
    [a, b, c] = deal(ELLIPSOIDE(4), ELLIPSOIDE(5), ELLIPSOIDE(6));

    normale = [(x - x0) / a^2; (y - y0) / b^2; (z - z0) / c^2];
    normale = normale / norm(normale);
end

function t_face = distanceIntersectionFace(pos, dir, plan, contraintes)
    [x0, y0, z0] = deal(pos(1), pos(2), pos(3));
    [vx, vy, vz] = deal(dir(1), dir(2), dir(3));
    [a, b, c, d] = deal(plan(1), plan(2), plan(3), plan(4));

    % Pour trouver l'intersection entre le rayon et le plan,
    % on plug l'equation parametrique du rayon dans l'equation du plan :
    % a*(x0 + t*vx) + b*(y0 + t*vy) + c*(z0 + t*vz) + d = 0
    % On obtient une equation du 1er degre en t
    A = a*vx + b*vy + c*vz;
    B = a*x0 + b*y0 + c*z0 + d;

    if A == 0 % Le rayon est parallele au plan
        t_face = Inf;
        return;
    end

    t = -B / A;
    if t <= 0 % Le plan est derriere le rayon
        t_face = Inf;
        return;
    end

    point_intersection = pos + t * dir;
    for i = 1:3
        limite_min = contraintes(i, 1);
        limite_max = contraintes(i, 2);
        if point_intersection(i) < limite_min || point_intersection(i) > limite_max
            t_face = Inf;
            return;
        end
    end

    t_face = t;
end

function R = coefficientReflexion(angle_incidence, angle_refraction, nin, nout)
    if nin == nout
        R_TE = 0;
        R_TM = 0;
    elseif angle_incidence == 0
        R_TE = ((nin - nout) / (nin + nout))^2;
        R_TM = R_TE;
    elseif abs(angle_incidence + angle_refraction) == pi/2
        R_TE = sin(angle_refraction - angle_incidence)^2;
        R_TM = 0;
    else
        R_TE = (sin(angle_refraction - angle_incidence) / sin(angle_refraction + angle_incidence))^2;
        R_TM = (tan(angle_refraction - angle_incidence) / tan(angle_refraction + angle_incidence))^2;
    end

    R = (R_TE + R_TM) / 2;
end