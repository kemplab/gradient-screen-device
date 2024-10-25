function Ed = MuSyC2(a, x1, x2)
% 2D Hill equation, drug properties and concentration

E1 = a(1); % Emax of drug 1
h1 = a(2); % Hill coefficient of drug 1
E0 = a(3); % E0 of drug 1 (and all)
C1 = a(4); % EC50 of drug 1

E2 = a(5); % Emax of drug 2
h2 = a(6); % Hill coefficient of drug 2
C2 = a(7); % EC50 of drug 2

alpha1 = a(8);
alpha2 = a(9);
% Assumes detailed balance

E3 = a(10); % beta is used in the initialization of MuSyC to calculate E3

r1 = a(11); % arbitrary number
r2 = a(12); % ^ not actually arbitrary in the full implementation
rn1 = (C1^h1)*r1; % this fulfills the condition C_i^(h_i) = (r_-i/r_i)
rn2 = (C2^h2)*r2;

conc1 = x1;
conc2 = x2;

% Make sure same size 
if size(conc1)==size(conc2)
else
    error('Concentration maps are not same size')
end

Ed = zeros(size(conc1));

for i = 1:numel(Ed) % loop through all elements sequentially

    d1 = conc1(i);
    d2 = conc2(i);

    % U
    Y11 = - r1.*d1.^h1 - r2.*d2.^h2;
    Y12 = rn1;
    Y13 = rn2;
    Y14 = 0;

    % A1
    Y21 = r1.*d1.^h1;
    Y22 = - rn1 - r2.*((alpha1.*d2).^h2);
    Y23 = 0;
    Y24 = rn2;

    % A2
    Y31 = r2.*d2.^h2;
    Y32 = 0;
    Y33 = - r1.*((alpha2.*d1).^h1) - rn2;
    Y34 = rn1;

    % A12
    Y41 = 1;
    Y42 = 1;
    Y43 = 1;
    Y44 = 1;

    Y1 = [Y11, Y12, Y13, Y14];
    Y2 = [Y21, Y22, Y23, Y24];
    Y3 = [Y31, Y32, Y33, Y34];
    Y4 = [Y41, Y42, Y43, Y44];

    Y = [Y1;Y2;Y3;Y4];

    Ev = [E0, E1, E2, E3];
    b = [0;0;0;1];
    Ed(i) = Ev/(Y)*b;
end

end