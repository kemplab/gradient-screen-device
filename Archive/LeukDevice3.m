%% Leukemia Drug Interactions Screening Device
% Non GUI version
% Normalized concentrations as multiples of EC50
%% Device Implementation
% Cell Seeding

tic;
conc1 = [0.5, 1, 1.5];
conc2 = [0.5, 1, 1.5];
conc3 = [0.5, 1, 1.5];
Emax1 = [0.25, 0.5];
hill1 = [1, 2, 4];
Emax2 = [0.25, 0.5];
hill2 = [1, 2, 4];
Emax3 = [0.25, 0.5];
hill3 = [1, 2, 4];
alpha12 = [-2, 0, 2];
alpha13 = [-2, 0, 2];
alpha23 = [-2, 0, 2];
beta12 = 0;
beta13 = 0;
beta23 = 0;

experimentDesign = fullfact([length(conc1), length(conc2), length(conc3), length(Emax1), length(hill1), ...
    length(Emax2), length(hill2), length(Emax3), length(hill3), ...
    length(alpha12), length(alpha13), length(alpha23), ...
    length(beta12), length(beta13), length(beta23)]);


experimentValues = experimentDesign;
for i = 1:height(experimentValues)
    experimentValues(i,1) = conc1(experimentDesign(i,1));
    experimentValues(i,2) = conc2(experimentDesign(i,2));
    experimentValues(i,3) = conc3(experimentDesign(i,3));
    experimentValues(i,4) = Emax1(experimentDesign(i,4));
    experimentValues(i,5) = hill1(experimentDesign(i,5));
    experimentValues(i,6) = Emax2(experimentDesign(i,6));
    experimentValues(i,7) = hill2(experimentDesign(i,7));
    experimentValues(i,8) = Emax3(experimentDesign(i,8));
    experimentValues(i,9) = hill3(experimentDesign(i,9));
    experimentValues(i,10) = alpha12(experimentDesign(i,10));
    experimentValues(i,11) = alpha13(experimentDesign(i,11));
    experimentValues(i,12) = alpha23(experimentDesign(i,12));
    experimentValues(i,13) = beta12(experimentDesign(i,13));
    experimentValues(i,14) = beta13(experimentDesign(i,14));
    experimentValues(i,15) = beta23(experimentDesign(i,15));
end

num = height(experimentDesign)
viability12 = zeros(num, 1);
viability13 = zeros(num, 1);
viability23 = zeros(num,1);


% Shut down
addpath('Concentration Maps')
data = load('ConcentrationMap_3uLpmin3D_D1e-10 - cropped.mat');
concMap = array2table(data.cc);
concMap.Properties.VariableNames = {'X', 'Y', 'c1', 'c2', 'c3'};
rmpath('Concentration Maps')
toc;

tic;
parfor i = 1:num
    channel_1_conc = experimentValues(i,1);
    channel_2_conc = experimentValues(i,2);
    channel_3_conc = experimentValues(i,3);
    Emax1 = experimentValues(i,4);
    hill1 = experimentValues(i,5);
    Emax2 = experimentValues(i,6);
    hill2 = experimentValues(i,7);
    Emax3 = experimentValues(i,8);
    hill3 = experimentValues(i,9);
    alpha12 = experimentValues(i,10);
    alpha13 = experimentValues(i,11);
    alpha23 = experimentValues(i,12);
    beta12 = experimentValues(i,13);
    beta13 = experimentValues(i,14);
    beta23 = experimentValues(i,15);
    
    [viability12(i), viability13(i), viability23(i)] = leukemiaScreen(concMap, ...
    channel_1_conc, channel_2_conc, channel_3_conc, ...
    Emax1, hill1, Emax2, hill2, Emax3, hill3, ...
    alpha12, alpha13, alpha23, beta12, beta13, beta23);
end
toc;


% Do studies looking only at viability12 compared to control drugs
synThreshold = 0.2; % 40%
synDiff = (min([viability23, viability13])-viability12)/(min([viability23, viability13]));
if synDiff > synThreshold
    synBool = true;
else
    synBool = false;
end



indvars = array2table(experimentValues, 'VariableNames',...
    {'conc1','conc2','conc3','Emax1','hill1','Emax2','hill2','Emax3','hill3',...
    'alpha12','alpha13','alpha23','beta12','beta13','beta23'});

results = [table(viability12), table(viability13), table(viability23), indvars];
writetable(results, 'results.csv');

function [viability12, viability13, viability23] = leukemiaScreen(concMap, ...
    channel_1_conc, channel_2_conc, channel_3_conc, Emax1, hill1, Emax2, hill2, Emax3, hill3, ...
    alpha12, alpha13, alpha23, beta12, beta13, beta23)

% Sample a subset of points
n = 3000; % number of cells
[cells, idx] = datasample(concMap, n);
%plot(cells.X, cells.Y,'.b')
%set(gca,'Color','k')
% Cell Properties

drug1_resist = 0; % percentage
drug2_resist = 0; 
drug3_resist = 0;
cell_death = 0;
% Channel Concentration of Drug
% Concentration relative to the EC50 (i.e. 1 = EC50, 2 = 2*EC50, etc.)

% Channel A is left of the top, B and C are clockwise away
%channel_1_conc = 1.6; % nM
%channel_2_conc = 1.4; % nM
%channel_3_conc = 0.8; % nM
cells.c1 = concMap.c1(idx).*(channel_1_conc);
cells.c2 = concMap.c2(idx).*(channel_2_conc);
cells.c3 = concMap.c3(idx).*(channel_3_conc);
%% Drug Properties

% Real drugs
% Single drug
%Emax1 = 0.35;
%hill1 = 1;
%Emax2 = 0.35;
%hill2 = 3;
%Emax3 = 0;
%hill3 = 3.6;
% Interactions
% Synergistic Potency

%alpha12 = 2;
%alpha23 = 1;
%alpha13 = 1;
%% 
% Synergistic Efficacy

%beta12 = 0;
%beta23 = 0;
%beta13 = 0;


drug1 = drugSingle(Emax1, hill1, 1, 1);
drug2 = drugSingle(Emax2, hill2, 1, 1);
drug3 = drugSingle(Emax3, hill3, 1, 1);
drug1.name = 'Drug 1';
drug2.name = 'Drug 2';
drug3.name = 'Drug 3';

% Pairwise combinations
drug12 = drugCombo(drug1, drug2);
drug12.alpha = [10^alpha12, 10^alpha12];
drug12.beta = beta12;

drug23 = drugCombo(drug2, drug3);
drug23.alpha = [10^alpha23, 10^alpha23];
drug23.beta = beta23;

drug13 = drugCombo(drug1, drug3);
drug13.alpha = [10^alpha13, 10^alpha13];
drug13.beta = beta13;

% Triple combination
drugCombination = drugTri(drug12, drug23, drug13);
drugCombination.gamma = [1, 1, 1];
drugCombination.drug1 = drug1;
drugCombination.drug2 = drug2;
drugCombination.drug3 = drug3;
%% Drug Effect

% Apply MuSyC model
Ed = MuSyC(drugCombination, cells.c1*(100-drug1_resist)/100, ...
    cells.c2*(100-drug2_resist)/100, cells.c3*(100-drug3_resist)/100);

% Kill cells according to Ed
live = rand(length(Ed),1)<Ed; % cells that are alive are 1, dead are 0

% Apply random cell death
rand_live = rand(length(Ed),1)>cell_death./100; % 1 means alive, 0 means dead
live = logical(live.*rand_live);

dead = live~=1; % reverse the 0 and 1 for live values
none = zeros(length(cells.c1),1); % blank matrix of image size
cells.live = live;
cells.dead = dead;

a_threshold = channel_1_conc.*0.1; % 10%
idx = find(cells.c1 < a_threshold);
cells_filtered23 = cells(idx,:);

b_threshold = channel_2_conc.*0.1; % 10%
idx = find(cells.c2 < b_threshold);
cells_filtered13 = cells(idx,:);

c_threshold = channel_3_conc.*0.1; % 10%
idx = find(cells.c3 < c_threshold);
cells_filtered12 = cells(idx,:);

viability12 = sum(cells_filtered12.live)/height(cells_filtered12);
viability13 = sum(cells_filtered13.live)/height(cells_filtered13);
viability23 = sum(cells_filtered23.live)/height(cells_filtered23);
viabilitycenter = cells.live;
viabilitycenter([idx_a; idx_b; idx_c]) = [];
viabilitycenter = sum(viabilitycenter)/height(viabilitycenter);

end