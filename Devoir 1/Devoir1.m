function [pcmNL, INL, alphaNL] = Devoir1(AngRot, vangulaire, forces, posNL)
    [pcmNL, INL, alphaNL] = NavetteLanceur(posNL, AngRot, vangulaire, forces);
end

function [pcmNL, INL, alphaNL] = NavetteLanceur(posNL, AngRot, vangulaire, forces)
    mN = 109e3; % Masse de la navette
    mR = 739e3; % Masse du réservoir
    mP = 469e3; % Masse du propulseur
    mNL = mN + mR + 2 * mP; % Masse du système

    [pcmN, IN] = Navette();
    [pcmR, IR] = Reservoir();
    [pcmPG, pcmPD, IPG, IPD] = Propulseurs();
    pcmNL = (mN * pcmN + mR * pcmR + mP * (pcmPG + pcmPD)) / mNL; % Position du centre de masse du système

    IN = IN + mN * T(pcmN, pcmNL);
    IR = IR + mR * T(pcmR, pcmNL);
    IPG = IPG + mP * T(pcmPG, pcmNL);
    IPD = IPD + mP * T(pcmPD, pcmNL);
    INL = IN + IR + IPG + IPD; % Matrice d'inertie du système

    rN = 3.5; % Rayon de la navette
    rR = 4.2; % Rayon du réservoir
    rP = 1.855; % Rayon des propulseurs
    forceN = [0; 0; forces(1)]; % Force de la navette
    forcePG = [0; 0; forces(2)]; % Force du propulseur gauche
    forcePD = [0; 0; forces(3)]; % Force du propulseur droit
    mfN = cross([0; 0; 0] - pcmNL, forceN); % Moment de force de la navette
    mfPG = cross([-rR - rP; rN + rR; 0] - pcmNL, forcePG); % Moment de force du propulseur gauche
    mfPD = cross([rR + rP; rN + rR; 0] - pcmNL, forcePD); % Moment de force du propulseur droit
    mfNL = mfN + mfPG + mfPD; % Moment de force total
    vangulaire = Rx(AngRot)' * vangulaire; % Rotation
    alphaNL = inv(INL) * (mfNL + cross(INL * vangulaire, vangulaire)); % cas 1 et cas 2

    pcmNL = Rx(AngRot) * pcmNL + posNL; % Rotation + Translation    cas 2 seulement
    INL = Rx(AngRot) * INL * Rx(AngRot)'; % Rotation
    alphaNL = Rx(AngRot) * alphaNL; % Rotation
end

function [pcmN, IN] = Navette()
    h1N = 27.93; % Hauteur du cylindre
    h2N = 9.31; % Hauteur du cône

    rN = 3.5; % Rayon de la navette
    V1N = pi * rN^2 * h1N; % Volume du cylindre
    V2N = pi * rN^2 * h2N / 3; % Volume du cône
    VN = V1N + V2N; % Volume de la navette

    mN = 109e3; % Masse de la navette
    m1N = mN * V1N / VN; % Masse du cylindre
    m2N = mN * V2N / VN; % Masse du cône

    pcm1N = [0; 0; h1N / 2]; % Position du centre de masse du cylindre
    pcm2N = [0; 0; h1N + h2N / 4]; % Position du centre de masse du cône
    pcmN = (m1N * pcm1N + m2N * pcm2N) / mN; % Position du centre de masse de la navette

    I1N = Icylindre(m1N, rN, h1N);
    I1N = I1N + m1N * T(pcm1N, pcmN);
    I2N = Icone(m2N, rN, h2N);
    I2N = I2N + m2N * T(pcm2N, pcmN);
    IN = I1N + I2N; % Matrice d'inertie de la navette
end

function [pcmR, IR] = Reservoir()
    hR = 39.1 + 7.8; % Hauteur du réservoir
    h1R = (2 / 3) * hR; % Hauteur du cylindre contenant l'hydrogène
    h2R = 39.1 - h1R; % Hauteur du cylindre contenant l'oxygène
    h3R = 7.8; % Hauteur du cône

    rR = 4.2; % Rayon du réservoir
    V2R = pi * rR^2 * h2R; % Volume du cylindre contenant l'oxygène
    V3R = pi * rR^2 * h3R / 3; % Volume du cône
    VO = V2R + V3R; % Volume de l'oxygène

    m1R = 108e3; % Masse de l'hydrogène
    mO = 631e3; % Masse de l'oxygène
    m2R = mO * V2R / VO; % Masse du cylindre contenant l'oxygène
    m3R = mO * V3R / VO; % Masse du cône

    rN = 3.5; % Rayon de la navette
    pcm1R = [0; rN + rR; h1R / 2]; % Position du centre de masse du cylindre contenant l'hydrogène
    pcm2R = [0; rN + rR; h1R + h2R / 2]; % Position du centre de masse du cylindre contenant l'oxygène
    pcm3R = [0; rN + rR; h1R + h2R + h3R / 4]; % Position du centre de masse du cône
    pcmR = (m1R * pcm1R + m2R * pcm2R + m3R * pcm3R) / (m1R + m2R + m3R); % Position du centre de masse du réservoir

    I1R = Icylindre(m1R, rR, h1R);
    I1R = I1R + m1R * T(pcm1R, pcmR);
    I2R = Icylindre(m2R, rR, h2R);
    I2R = I2R + m2R * T(pcm2R, pcmR);
    I3R = Icone(m3R, rR, h3R);
    I3R = I3R + m3R * T(pcm3R, pcmR);
    IR = I1R + I2R + I3R; % Matrice d'inertie du réservoir
end

function [pcmPG, pcmPD, IPG, IPD] = Propulseurs()
    h1P = 39.9; % Hauteur des cylindres
    h2P = 5.6; % Hauteur des cônes

    rP = 1.855; % Rayon des propulseurs
    V1P = pi * rP^2 * h1P; % Volume du cylindre
    V2P = pi * rP^2 * h2P / 3; % Volume du cône
    VP = V1P + V2P; % Volume de chaque propulseur

    mP = 469e3; % Masse du propulseur
    m1P = mP * V1P / VP; % Masse du cylindre
    m2P = mP * V2P / VP; % Masse du cône

    rN = 3.5; % Rayon de la navette
    rR = 4.2; % Rayon du réservoir
    pcm1PG = [-rR - rP; rN + rR; h1P / 2]; % Position du centre de masse du cylindre du propulseur gauche
    pcm1PD = [rR + rP; rN + rR; h1P / 2]; % Position du centre de masse du cylindre du propulseur droit
    pcm2PG = [-rR - rP; rN + rR; h1P + h2P / 4]; % Position du centre de masse du cône du propulseur gauche
    pcm2PD = [rR + rP; rN + rR; h1P + h2P / 4]; % Position du centre de masse du cône du propulseur droit
    pcmPG = (m1P * pcm1PG + m2P * pcm2PG) / mP; % Position du centre de masse du propulseur gauche
    pcmPD = (m1P * pcm1PD + m2P * pcm2PD) / mP; % Position du centre de masse du propulseur droit

    I1PG = Icylindre(m1P, rP, h1P);
    I1PG = I1PG + m1P * T(pcm1PG, pcmPG);
    I2PG = Icone(m2P, rP, h2P);
    I2PG = I2PG + m2P * T(pcm2PG, pcmPG);
    I1PD = Icylindre(m1P, rP, h1P);
    I1PD = I1PD + m1P * T(pcm1PD, pcmPD);
    I2PD = Icone(m2P, rP, h2P);
    I2PD = I2PD + m2P * T(pcm2PD, pcmPD);
    IPG = I1PG + I2PG; % Matrice d'inertie du propulseur gauche
    IPD = I1PD + I2PD; % Matrice d'inertie du propulseur droit
end

function Rx = Rx(AngRot)
    Rx = [
        1, 0, 0;
        0, cos(AngRot), -sin(AngRot);
        0, sin(AngRot), cos(AngRot)
    ];
end

function T = T(pcm, p)
    t = p - pcm;
    Txx = t(2)^2 + t(3)^2;
    Txy = -t(1) * t(2);
    Txz = -t(1) * t(3);
    Tyx = -t(2) * t(1);
    Tyy = t(1)^2 + t(3)^2;
    Tyz = -t(2) * t(3);
    Tzx = -t(3) * t(1);
    Tzy = -t(3) * t(2);
    Tzz = t(1)^2 + t(2)^2;
    T = [Txx, Txy, Txz; Tyx, Tyy, Tyz; Tzx, Tzy, Tzz];
end

function Icylindre = Icylindre(m, r, h)
    Ixx = (m / 4) * (r^2) + (m / 12) * (h^2); % Moment d'inertie du cylindre par rapport à l'axe x
    Iyy = Ixx; % Moment d'inertie du cylindre par rapport à l'axe y
    Izz = (m / 2) * (r^2); % Moment d'inertie du cylindre par rapport à l'axe z
    Icylindre = [Ixx, 0, 0; 0, Iyy, 0; 0, 0, Izz]; % Matrice d'inertie du cylindre
end

function Icone = Icone(m, r, h)
    Ixx = (m / 80) * (12 * r^2 + 3 * h^2); % Moment d'inertie du cône par rapport à l'axe x
    Iyy = Ixx; % Moment d'inertie du cône par rapport à l'axe y
    Izz = (m / 10) * (3 * r^2); % Moment d'inertie du cône par rapport à l'axe z
    Icone = [Ixx, 0, 0; 0, Iyy, 0; 0, 0, Izz]; % Matrice d'inertie du cône
end