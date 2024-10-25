clear;clc;close all
addpath(['..', filesep, 'DataDerived'])


T1 = readtable('20210920_dev1.csv');
T2 = readtable('20210920_dev2.csv');
T3 = readtable('20210920_dev3.csv');

T = [T1; T2; T3]; % tables all stacked?

% rename variables to match cells format
T = renamevars(T,["x","y","C1","C2","C3"], ...
                 ["X","Y","c1","c2","c3"]);
  

% load ternary plot packages
addpath('TernaryPlot')
addpath('ternary2')


% create zones
cells = T;
thresholdMainZone = 0.08; % below __% of drug is considered no drug
cells.zone = zeros(height(cells), 1);
cells.viability = zeros(height(cells), 1);

% zone 1: interaction between drug 1 and 2 (no 3)
% get threshold based on max
idx{1} = find(cells.c1 < thresholdMainZone);
idx{2} = find(cells.c2 < thresholdMainZone);
idx{3} = find(cells.c3 < thresholdMainZone);

% assign zone number in table
for i = 1:length(idx)
    cells.zone(idx{i}) = i;
end
% zone 0: middle portion of device
% zone 1: interaction between 2 and 3
% zone 2: interaction between 1 and 3
% zone 3: interaction between 1 and 2

% perform a moving window average across data (using convolution function conv)
% test examples a and b (for no drug c)

% viability along drug23
a = cells.c2(idx{1}); % concentration of a
b = cells.c3(idx{1}); % concentration of b
live = cells.live(idx{1});
sumab = cells.c2(idx{1}) + cells.c3(idx{1}); % check to see total proportion of drug
xaxis = b./sumab; % normalize b to the sumab


filter_width = 0.25;

% create new xaxis [0, 1] proportion of b to sumab
[sortedAxis, I] = sort(xaxis);
liveSorted = live(I);
viability = movmean(liveSorted, length(idx{1}).* filter_width); % 20% of windows sorting

figure
plot(sortedAxis, viability, '-') % moving average across __% of cells
[revSort, I2] = sort(I); % reverse sort
cells.viability(idx{1}) = viability(I2); % assign viability to main table

x4fit_d1 = a./max(a); % normalize
y4fit_d2 = b./max(b);
z = viability(I2);


ED50_1 = 2.7;
ED50_2 = 2.7;
E0 = 1;
Emax = 0;
h1 = 2.5;
h2 = 2.5;

syms E d1 d2 a
termA = (d1 ./ (ED50_1 .* ((E - E0) ./ (Emax - E)).^ (1./h1)));
termB = (d2 ./ (ED50_2 .* ((E - E0) ./ (Emax - E)).^ (1./h2)));
termC = a .* d1 .* d2;
termD = ED50_1 .* ED50_2 .* ((E - E0) ./ (Emax - E)) .^ ((1 ./ (2 .* h1)) + (1 ./ (2 .* h2)));
equation = termA + termB + (termC ./ termD); % -1 to solve for equal to 0

test = solve(equation == 1, E, 'Real', true, 'IgnoreAnalyticConstraints', true);

eqn2solve = char(test);

ft = fittype( eqn2solve, ...
    'independent', {'d1', 'd2'}, 'dependent', 'E')

d1_max = 1;
d2_max = 1;
d1 = linspace(0.001, d1_max, 100); % left side
d2 = linspace(d2_max, 0.001, 100); % right side
xaxis4fit = linspace(0, 1, 100)

opts = fitoptions( ft );
opts.Lower = [-10];
%opts.Upper = [10];
opts.Robust = 'LAR';
opts.StartPoint = [1];

[f, gof] = fit( [x4fit_d1, y4fit_d2], z, ft,...
    opts)

% plot the fit on top
a = f.a;
y_fitted = eval(subs(test));
hold on
plot(xaxis4fit, y_fitted,'--')
hold off

figure
plot(f, [x4fit_d1, y4fit_d2], z)