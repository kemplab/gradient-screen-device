%% Plot CellProfiler points on rotated and scaled concentration map
% Generate concentration map
grid100 = csvread('5umgrid.csv',1,0);
concentrationMap = generateConcentrationMap(grid100);
mapcentroid = concentrationMapCentroid(concentrationMap);

% Import all original images
imageraw = imread('CalceinEthd1Hoechst_Overlay.tif');

mapscale = sum(round(sum(concentrationMap(1,:,:),3)));
channelWidth = 540.0316; % from new method
scalefactor = channelWidth/mapscale;

% Set rotation and center variables
%imangle = theta*180/pi; % degrees, CCW, this is coded to be 39.6606 deg CCW
imangle = 121.5889-90;
imscale = scalefactor;
imcenter = mapcentroid;

% EXECUTE CONCENTRATION MAP
	concentrationMapRSC = RotateScaleCrop(concentrationMap, imscale, imangle, imcenter);

% Hard coded image center
imagecenter = round(intersection); % pixel needs to be integer
% making Y of center smaller, moves the map DOWN
falseimage = zeros(3066,4086);
croppedMap = cropMap2Image(concentrationMapRSC, imageraw, imagecenter);


% Import cell data
calceinData = readtable('20200822_Live_Secondary.csv');
histogram(calceinData.Intensity_MeanIntensity_CropLive)
thresh = 0.02;

drug1C = croppedMap(:,:,1);
% find x and y locations of drug concentrations on image
x = calceinData.Location_Center_X;
y = calceinData.Location_Center_Y;

x = round(x);
y = round(y);

x = x + 1020; % to account for cellProfiler cropping
y = y + 240;
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

% figure(6) plot live/dead simulated
figure(6); clf
imshow(falseimage)
hold on
for i = 1:length(live)
    if live(i)==1
        plot(x(i),y(i),'.g')
    else
        plot(x(i),y(i),'.r')
    end
end
hold off

figure(7)
%semilogx(C1(find(live==1)), viability(find(live==1)), '.')
semilogx(C1, viability, '.')
xlabel('Concentration')
ylabel('Viability')