% Creates "smile plots" for evaluation of pairwise drug synergy
% from continuous screen data from microfluidic device
% Author: Dan Zhang
% Updated: 2021-12-16


clear;clc;close all

% from microfluidic device experiments, each dev is a replicate

% Drug ordering
% Drug 1: Methotrexate
% Drug 2: MRX-2843
% Drug 3: Vincristine

% colors for drug pairs (RGB, normalized to matlab [0 1] scaling)
% Ordering of pairs: 2-3, 1-3, 1-2
color.darkBlue = [0, 63, 92]/255;
color.magenta = [188, 80, 144]/255;
color.orange = [255, 166, 0]/255;
color.grey = [100, 100, 100]/255;


T = readtable('20210422_Dev12_concentrations.csv');

% load ternary plot packages
addpath('TernaryPlot')
addpath('ternary2')

% create zones
cells = T;
thresholdMainZone = 0.08; % below 8% of drug is ignored
cells.zone = zeros(height(cells), 1);
cells.viability = zeros(height(cells), 1);
[cells.ternX, cells.ternY] = ternCoord(cells.c1, cells.c2, cells.c3);

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


%------------------------------------------------------------------------------
% Figure 1. Pairwise lines plotted against ratio B/total
figure
hold on
filter_width = 0.2;
[T23, T13, T12] = isolatePair(T, 8); % 8% of drug
[x23, viability23] = viabilityOverDrugRatio(T23, filter_width);
[x13, viability13] = viabilityOverDrugRatio(T13, filter_width);
[x12, viability12] = viabilityOverDrugRatio(T12, filter_width);

viability23_smoothed = medfilt1(smooth(viability23, .2));
viability13_smoothed = medfilt1(smooth(viability13, .2));
viability12_smoothed = medfilt1(smooth(viability12, .2));

plot(x23, viability23_smoothed, '-', 'LineWidth', 3, 'Color', color.darkBlue)
plot(x13, viability13_smoothed, '-', 'LineWidth', 3, 'Color', color.magenta)
plot(x12, viability12_smoothed, '-', 'LineWidth', 3, 'Color', color.orange)


% Lines of additivity
% Plot from first point to last point
%plot([x23(1), x23(end)], [viability23(1), viability23(end)], '-', 'Color', color.grey, 'LineWidth', 2)
%plot([x13(1), x13(end)], [viability13(1), viability13(end)], '-', 'Color', color.grey, 'LineWidth', 2)
%plot([x12(1), x12(end)], [viability12(1), viability12(end)], '-', 'Color', color.grey, 'LineWidth', 2)

% legend and chart settings
legend({"MRX-VCR", "MTX-VCR", "MTX-MRX", "Lines of Additivity"}, 'Location', 'SouthEast')
xlabel('Ratio Drug B/A')
ylabel('Viability')
set(gca,'FontSize',14)
axis([0 1 0 1])

%% CUSTOM FUNCTIONS ================================================================

% -----------------------------------------------------------------------------
function AUC = viabilityAUC(T)
% Calculate the AUC of the viability data, assuming straight line from left and right sides
% some y=mx+b stuff happening

% create new axis
concentrations = T{:, 3:5};
if mean(concentrations(:,1)) < 0.1
    a = concentrations(:,2);
    b = concentrations(:,3);
elseif mean(concentrations(:,2)) < 0.1
    a = concentrations(:,1);
    b = concentrations(:,3);
elseif mean(concentrations(:,3)) < 0.1
    a = concentrations(:,1);
    b = concentrations(:,2);
end

xaxis = b./(a + b);
[sortedAxis, I] = sort(xaxis);
liveSorted = T.live(I);
filter_width = 0.1;
viability = movmean(liveSorted, height(T).* filter_width);
viability = medfilt1(smooth(viability, .2)); % doubly smoothed line
x = sortedAxis;
y = viability;

% finding the line of additivity
x_line = linspace(x(1), x(end), 20); % creating new x for the straight line
m = (y(end) - y(1))./(x(end) - x(1));
b = y(1) - m.*(x(1));
y_line = m.*x_line + b;

% subtracting AUCs
AUC_additivity = trapz(x_line, y_line);
AUC_data = trapz(x, y);
AUC = AUC_additivity - AUC_data;

% Plotting script
%{
figure
plot(x, y, '-')
hold on
plot(x_line, y_line, '.--')
hold off
axis([0, 1, 0, 1])
%}
end % end of function

% -----------------------------------------------------------------------------
function viability = avgViability(T)
viability = mean(T.live); 
end % end of function

% -----------------------------------------------------------------------------
% isolatePair: extract pairwise data under cutoff
function [T23, T13, T12] = isolatePair(T, cutoffPercent)
% taking a table of cell data, produce only the data with under the percent cutoff (%)
thresholdMainZone = cutoffPercent./100;
cells = T;
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
T23 = cells(idx{1}, :);
T13 = cells(idx{2}, :);
T12 = cells(idx{3}, :);
end % end of function

% -----------------------------------------------------------------------------
% viabilityOverDrugRatio: creates normalized ratiometric axis by viability
function [x, viability] = viabilityOverDrugRatio(T, windowWidth)
% windowWidth defines moving average filter width, max 1

% which drug is missing? 1, 2, or 3
[~, missing] = min([mean(T.c1), mean(T.c2), mean(T.c3)]);
switch missing
case 1
    druga = T.c2;
    drugb = T.c3;
case 2
    druga = T.c1;
    drugb = T.c3;
case 3
    druga = T.c1;
    drugb = T.c2;
end

% create new ratiometric x axis
live = T.live;
sumab = druga + drugb;
xaxis = drugb./sumab;

% sort the axis for averaging
[x, I] = sort(xaxis);
liveSorted = live(I);
viability = movmean(liveSorted, height(T)*windowWidth);
end % end of function

% -----------------------------------------------------------------------------
% ternCoord: convert cartesian to ternary coordinate system
function [ternX, ternY] = ternCoord(a,b,c)
    ternX = 0.5.*(2.*b+c)./(a+b+c);
    ternY = (sqrt(3)./2).*c./(a+b+c);
end % end of function

% -----------------------------------------------------------------------------
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
end % end of function

% -----------------------------------------------------------------------------
function viabilityOriginalSort = viabilityPairTernary(conc1, conc2, live)
% for a linear combination of concentrations of drugs 1 and 2, calculates viability
% as a moving window with increasing ratio of drug 1/2
sum12 = conc1 + conc2;
xaxis = conc2./sum12; % normalize concentration of drug 2
% (so that when sorted, goes from drug 1 to drug 2)

% sort
[sortedAxis, I] = sort(xaxis);
liveSorted = live(I);
viability = movmean(liveSorted, length(conc1).*0.20); % moving average of 20% of cells

[sorted2 sortIndRev] = sort(I);
viabilityOriginalSort = viability(sortIndRev);
end % end of function

% -----------------------------------------------------------------------------
function binnedAUC = ratioSpecificAUC(T)

% create new axis

% sort by new axis


% bin ratios (by 10%) and live dead


% viability from additivity assumption


% ratiometric synergy calculation

% plot bar plot of bins

end % end of function