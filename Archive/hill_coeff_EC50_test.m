% question: is a change in hill same for differing levels of EC50?
drug1 = struct;
drug2 = struct;

E0 = 1;
Emax = 0.2;
EC50_1 = 0.0050;
EC50_2 = 0.1;
hill_coeff = 8;

% set drug 1 parameters
drug1.E0 = E0;
drug1.Emax = Emax;
drug1.C = EC50_1;
drug1.h = hill_coeff;

% set drug 2 parameters
drug2.E0 = E0;
drug2.Emax = Emax;
drug2.C = EC50_2;
drug2.h = hill_coeff;

d1 = linspace(0, 0.05, 11); % set original dose points
d2 = linspace(0, 1, 11);

% run test
E1 = hill(d1, drug1);
E2 = hill(d2, drug2);

% conclusion: hill coefficient is just a shape parameter so everything
% will look the same normalized to a consistent EC50 dosing

function E = hill(d, params)
    E0 = params.E0;
    Emax = params.Emax;
    C = params.C;
    h = params.h;
    
    E = E0 + (Emax - E0) .* (d.^h ./ (C.^h + d.^h));
end