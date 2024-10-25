function Ed = MuSyC1(a,x)
% 1D Hill equation, coefficients and concentration
% a is coefficients vector
% x is independent variable (drug concentration)

Emax = a(1); % Maximum effect of the drug
h = a(2); % Hill coefficient
E0 = a(3); % Effect at 0 drug
C = a(4); % EC50

Ed = ((E0-Emax).*C.^h)./(C.^h + x.^h) + Emax;
end