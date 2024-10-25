%% set drug parameters
% Single drug
Emax1 = 0;
hill1 = 3;
Emax2 = 0;
hill2 = 3;
Emax3 = 0;
hill3 = 3;

% Interactions
% Synergistic Potency
alpha12 = 0;
alpha23 = 0;
alpha13 = 0;

% Synergistic Efficacy
beta12 = 0;
beta23 = 0;
beta13 = 0;

drug1 = DrugSingle(Emax1, hill1, 1, 1);
drug2 = DrugSingle(Emax2, hill2, 1, 1);
drug3 = DrugSingle(Emax3, hill3, 1, 1);
drug1.Name = 'Drug 1';
drug2.Name = 'Drug 2';
drug3.Name = 'Drug 3';

% Pairwise combinations
drug12 = DrugCombo(drug1, drug2);
drug12.Alpha = [10^alpha12, 10^alpha12];
drug12.Beta = beta12;

drug23 = DrugCombo(drug2, drug3);
drug23.Alpha = [10^alpha23, 10^alpha23];
drug23.Beta = beta23;

drug13 = DrugCombo(drug1, drug3);
drug13.Alpha = [10^alpha13, 10^alpha13];
drug13.Beta = beta23;

% Triple combination
drugCombination = DrugTri(drug12, drug13, drug23);
drugCombination.Gamma = [1, 1, 1];
drugCombination.Drug1 = drug1;
drugCombination.Drug2 = drug2;
drugCombination.Drug3 = drug3;

%% Defining input concentration data
c1 = 1;
c2 = 1;
c3 = 1;

%% Applying MuSyC model
Ed = MuSyC(drugCombination, c1, ...
    c2, c3)

