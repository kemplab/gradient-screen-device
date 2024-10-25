
%{
syms IC50_1 IC50_2 drug_1 drug_2 E E_CON a m_1 m_2
termA = (drug_1 ./ (IC50_1 .* (E ./ (E_CON - E)) .^ (1 ./ m_1)));
termB = (drug_2 ./ (IC50_2 .* (E ./ (E_CON - E)) .^ (1 ./ m_2)));
termC = a .* drug_1 .* drug_2;
termD = IC50_1 .* IC50_2 .* (E ./ (E_CON - E)) .^ ((1 ./ (2 .* m_1)) + (1 ./ (2 .* m_2)));
equation = termA + termB + termC ./ termD == 1; % -1 to solve for equal to 0



% Fill in simulated parameters to calculate a
IC50_1 = 1;
IC50_2 = 1;
m_1 = 2;
m_2 = 1;

E_CON = 1;
a = 4;
% assumed direction for synergy
% a > 0 = synergy
% a < 0 = antagonism

drug_1 = 0.3;
drug_2 = 0.3;


eq2 = solve(subs(equation), E, 'Real', true);
eval(eq2)


drug_1 = [0.001 1 2 3 4 5];
test2 = solve(termA==1, E, 'Real', true)

y = eval(subs(test2));
plot(drug_1, y, 'o-')

%}


%{
% Alternative formulation model form of G L Drusano et al
syms d1 d2 ED50_1 ED50_2 E E0 Emax h1 h2 a
termA = (d1 ./ (ED50_1 .* ((E - E0) ./ (Emax - E)).^ (1./h1)));
termB = (d1 ./ (ED50_2 .* ((E - E0) ./ (Emax - E)).^ (1./h2)));
termC = a .* d1 .* d2;
termD = ED50_1 .* ED50_2 .* ((E - E0) ./ (Emax - E)) .^ ((1 ./ (2 .* h1)) + (1 ./ (2 .* h2)));
equation = termA + termB + (termC ./ termD); % -1 to solve for equal to 0


ED50_1 = 1;
ED50_2 = 1;
E0 = 1;
Emax = 0;
h1 = 3;
h2 = 3;
d1 = linspace(0.001, 5, 100);
d2 = linspace(5, 0.001, 100);
a = 1;

test3 = vpasolve(termA == 1, E, [0 1]);
% test3 = solve(termA == 1, E, 'Real', true); % using the normal solve
y = eval(subs(test3));
figure
plot(d1, y, '.-')

%test4 = vpasolve(equation == 1, E);
%y2 = eval(subs(test4));
%}



%% reducing the number of symbolic commands in the expression
d1 = linspace(0.001, 3, 20);
d2 = linspace(0.001, 3, 20);
y = test(d1, d2, 1);
plot(d1, y, '.-')
surf(real(y))
view(0, 90) % top down view

d1 = 1;
d2 = 1;
y = test(d1, d2, 1)


% actual function space

function y = test(d_1, d_2, a)
% a is interaction index from Greco et al. 1995
for i = 1:length(d_1);
    for j = 1:length(d_2);
ED50_1 = 1;
ED50_2 = 1;
E0 = 1;
Emax = 0;
h1 = 3;
h2 = 3;
d1 = d_1(i);
d2 = d_2(j);

% greco 1995 model, formulation of G L Drusano et al
syms E
termA = (d1 ./ (ED50_1 .* ((E - E0) ./ (Emax - E)).^ (1./h1)));
termB = (d2 ./ (ED50_2 .* ((E - E0) ./ (Emax - E)).^ (1./h2)));
termC = a .* d1 .* d2;
termD = ED50_1 .* ED50_2 .* ((E - E0) ./ (Emax - E)) .^ ((1 ./ (2 .* h1)) + (1 ./ (2 .* h2)));
equation = termA + termB + (termC ./ termD); % -1 to solve for equal to 0

test = vpasolve(equation == 1, E);
y(i,j) = eval(subs(test));

end % end of j loop
end % end of i loop
end % end of function
%test4 = vpasolve(equation == 1, E);
%y2 = eval(subs(test4));


% notes
% currently works for the vpasolve implementation of just d1
% trying to do this with the implementation for the entire equation