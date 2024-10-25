clear;clc; close all

% Input settings
selectFile = 'Dev11';

switch selectFile
    case 'Dev11'
        intersection = [1726-10; 1418-90];
        angleRotate = 31;
        meanChannelWidth = 392;
        filename = ['20210422_', selectFile];
    case 'Dev12'
        intersection = [1726+50; 1418-50];
        angleRotate = 34;
        meanChannelWidth = 392;
        filename = ['20210422_', selectFile];
    case 'Dev13'
        intersection = [1726+70; 1418-50];
        angleRotate = 29;
        meanChannelWidth = 392;
        filename = ['20210422_', selectFile];
    case 'Dev14'
        intersection = [1726+10; 1418-0];
        angleRotate = 29;
        meanChannelWidth = 392;
        filename = ['20210422_', selectFile];
    case 'Dev16'
        intersection = [1726+45; 1418-60];
        angleRotate = 33;
        meanChannelWidth = 392;
        filename = ['20210422_', selectFile];
    case 'Dev17'
        intersection = [1726-10; 1418-25];
        angleRotate = 31;
        meanChannelWidth = 392;
        filename = ['20210422_', selectFile];
end


%% Plot CellProfiler points on rotated and scaled concentration map
% Generate concentration map
grid100 = csvread('xy_25um_c1c2c3.csv');
gridcenter = [mean(grid100(:,1)), mean(grid100(:,2))]; % centroid as center
coord = grid100(:,1:2); % x and y coordinates

% scale
coord(:,1) = coord(:,1) - gridcenter(1); % normalize x
coord(:,2) = coord(:,2) - gridcenter(2); % normalize y
x = grid100(:,1);
y = grid100(:,2);
A = x(find(y==24675));
mapscale = max(A) - min(A);
scalefactor = meanChannelWidth/mapscale;
coord = coord.*scalefactor;

% rotate
theta = angleRotate; % theta rotates counterclockwise, degrees
R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
coord = coord*R';

% adjust new center
coord(:,1) = coord(:,1) + intersection(1); % add image center
coord(:,2) = coord(:,2) + intersection(2); % add image center


figure
plot(coord(:,1), coord(:,2),'.')


% Import cell data
calceinData = readtable([filename, 'Calcein', '.csv']);
PI = readtable([filename, 'PI', '.csv']);


%drug1C = croppedMap(:,:,1);
% find x and y locations of drug concentrations on image
x = calceinData.Location_Center_X;
y = calceinData.Location_Center_Y;
hold on
plot(x,y,'.r')

points = [x, y];
k = dsearchn(coord,points);
conc = grid100(k,3:5);

c1 = conc(:,1);
c2 = conc(:,2);
c3 = conc(:,3);
calceinMean = calceinData{:,11};
calceinMedian = calceinData{:,12};
PIMean = PI{:,11};
PIMedian = PI{:,12};
results = table(x,y,c1,c2,c3,calceinMean,calceinMedian,PIMean,PIMedian);
writetable(results,[filename, '_concentrations.csv']);
saveas(gcf, [filename, '_plot.png']);
%{
thresh = 0.01;

calceinValue = calceinData.Intensity_MeanIntensity_CropLive;
live = calceinValue>thresh;


% Rolling ball method for viability based on cells within radius r
% For each x,y cell, find indices of closest x,y within r
viability = x; % initilize variable
r = 500; % range for rolling
Idx = rangesearch([x,y],[x,y],r);
for i = 1:length(Idx) % for each cell
	% average live status of neighbors
	viability(i) = mean(live(Idx{i}));
end

C1 = x;
C2 = C1;
C3 = C1;
for i = 1:length(x)
    C1(i) = croppedMap(y(i),x(i),1); % index must be flipped
    C2(i) = croppedMap(y(i),x(i),2);
    C3(i) = croppedMap(y(i),x(i),3);
end

results = table(x,y,C1,C2,C3,viability,live,calceinValue);
%% FIGURES

% Figure 2. Default concentration map.
figure(2); clf
imagesc(concentrationMap(:,:,1))
hold on
plot(mapcentroid(1), mapcentroid(2),'r+')
hold off
title('Concentration Map with Center Point')


% Figure 3. Rotated, scaled, and centered, map.
figure(3); clf
imagesc(concentrationMapRSC(:,:,1))

% Figure 4. Cropped image
figure(4); clf
imagesc(croppedMap(:,:,1))

hold on
for i = 1:length(live)
    if live(i)==1
        plot(x(i),y(i),'.k')
    else
        plot(x(i),y(i),'.r')
    end
end
hold off


% Figure 5. Image of calcein with points highlighted
figure(5); clf
imshow(imageraw)
hold on
%for i = 1:cell_N
%	plot(cell_location(i,1), cell_location(i,2),'.g')
%end
hold off

%}