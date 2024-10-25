% Divides cell data into ternary zones and plots subzones
% Each zone is tagged with a number for subanalysis
function cells = ternaryZoning(cells)
% Cells structure:
% cols: X, Y, c1, c2, c3, live, dead

% load ternary plot packages
addpath('TernaryPlot')
addpath('ternary2')


% create zones

thresholdMainZone = 0.08; % below __% of drug is considered no drug
cells.zone = zeros(height(cells), 1);
cells.viability = zeros(height(cells), 1);
[cells.ternX cells.ternY] = ternCoord(cells.c1, cells.c2, cells.c3);

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

figure
tersurf(cells.c1, cells.c2, cells.c3, cells.live)
%==============================================================================
% figure 1. zones highlighted on ternary plot
f = figure;
tiledlayout(1,2)

% disable plotting on spatial grid
%{
ax1 = nexttile;
plot(cells.X, cells.Y, '.k')
hold on
for i = 1:length(idx)
    plot(cells.X(idx{i}), cells.Y(idx{i}), '.')
end
%}
ax2 = nexttile;
plot(cells.ternX, cells.ternY, '.k')
hold on
for i = 1:length(idx)
    plot(cells.ternX(idx{i}), cells.ternY(idx{i}), '.')
end
axis([0, 1, 0, 0.5*tand(60)])

f.Position = [300 300 1200 420];



%==============================================================================
% figure 2. "smile plot" viability axis chart thing

% perform a moving window average across data (using convolution function conv)
% test examples a and b (for no drug c)

% viability along drug23
a = cells.c2(idx{1}); % concentration of a
b = cells.c3(idx{1}); % concentration of b
live = cells.live(idx{1});
sumab = cells.c2(idx{1}) + cells.c3(idx{1}); % check to see total proportion of drug
xaxis = b./sumab; % normalize b to the sumab

% create new xaxis [0, 1] proportion of b to sumab
[sortedAxis, I] = sort(xaxis);
liveSorted = live(I);
viability = movmean(liveSorted, length(idx{1}).*0.20); % 20% of windows sorting

figure
plot(sortedAxis, viability, '-') % moving average across __% of cells
[revSort, I2] = sort(I); % reverse sort
cells.viability(idx{1}) = viability(I2); % assign viability to main table

% viability along drug13
a = cells.c1(idx{2});
b = cells.c3(idx{2});
live = cells.live(idx{2});
sumab = cells.c1(idx{2}) + cells.c3(idx{2}); % check to see total proportion of drug
xaxis = b./sumab;

[sortedAxis, I] = sort(xaxis);
liveSorted = live(I);
viability = movmean(liveSorted, length(idx{1}).*0.20);
hold on
plot(sortedAxis, viability, '-') % moving average across 10% of cells
[revSort, I2] = sort(I); % reverse sort
cells.viability(idx{2}) = viability(I2); % assign viability to main table


% viability along drug12
a = cells.c1(idx{3});
b = cells.c2(idx{3});
live = cells.live(idx{3});
sumab = cells.c1(idx{3}) + cells.c2(idx{3}); % check to see total proportion of drug
xaxis = b./sumab;

[sortedAxis, I] = sort(xaxis);
liveSorted = live(I);
viability = movmean(liveSorted, length(idx{1}).*0.20);
plot(sortedAxis, viability, '-') % moving average 20% of cells meeting criteria
axis([0 1 0 1])
[revSort, I2] = sort(I); % reverse sort
cells.viability(idx{3}) = viability(I2); % assign viability to main table
title('Viability along drugAB axis')

legend({"Drug23 axis", "Drug13 axis", "Drug12 axis"}, 'Location', 'SouthEast')
xlabel('Proportion of drug A to B')
ylabel('Viability')


% Rolling ball for center viability
idx{4} = find(cells.zone == 0);
%cells.viability(idx{4}) = rollingViability(cells(idx{4},:), 400); % what units is this r?
% ^ temporariliy disabled, rolling ball viability


%==============================================================================
% figure 3. plot ternary plot with shaded zones averaged across
figure
tersurf(cells.c1, cells.c2, cells.c3, cells.viability)

% ignore rolling ball viability
%{
figure % only rolling ball, compare accuracy
cells2 = cells;
cells2.viability = rollingViability(cells2, 400);
tersurf(cells2.c1, cells2.c2, cells2.c3, cells2.viability)
%}



%% CUSTOM FUNCTIONS ================================================================
function [ternX, ternY] = ternCoord(a,b,c)
    ternX = 0.5.*(2.*b+c)./(a+b+c);
    ternY = (sqrt(3)./2).*c./(a+b+c);
end


function viability = rollingViability(cells, r)
    % Rolling ball method for viability based on cells within radius r
    x = cells.X;
    y = cells.Y;
    live = cells.live;
    viability = cells.X; % initialize
    
    % For each x,y cell, find indices of closest x,y within r
    Idx = rangesearch([x,y],[x,y],r);
    
    for i = 1:length(Idx) % for each cell
        % average live status of neighbors
        viability(i) = mean(live(Idx{i}));
    end    
end


function viabilityOriginalSort = viabilityPairTernary(conc1, conc2, live)
% for a linear combination of concentrations of drugs 1 and 2, calculates viability
% as a moving window with increasing ratio of drug 1/2
sum12 = conc1 + conc2;
xaxis = conc2./sum12; % normalize concentration of drug 2 (so that when sorted, goes from drug 1 to drug 2)

% sort
[sortedAxis, I] = sort(xaxis);
liveSorted = live(I);
viability = movmean(liveSorted, length(conc1).*0.20); % moving average of 20% of total cells

[sorted2 sortIndRev] = sort(I);
viabilityOriginalSort = viability(sortIndRev);
end




end