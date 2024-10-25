
% Generate concentration map
grid100 = csvread('5umgrid.csv',1,0);
concentrationMap = generateConcentrationMap(grid100);
mapcentroid = concentrationMapCentroid(concentrationMap);


% Import all original images
%annexin = imread('annexin.tif'); % calcein
%calcein = imread('calcein.tif'); % annexin
imageraw = imread('CalceinEthd1Hoechst_Overlay.tif');


% hard coded angle and scaling from image
%line1 = [931, 1412; 2312, 3045];
%line2 = [1916, 718; 3073, 2141];
%normdist = lineNormDist(line1, line2);

% unsure if this part is needed but getting anyway
%lines = [1.64953271028037,-1362.83271759606;1.60283687943262,-2326.20690205549;-0.0192616372391654,1333.74866107748;-0.0312500000000000,1917.01541095890;-1.82795698924731,4857.60455989731;-1.76751592356688,5877.69093322697];



%perpline = [1727, 2371; 2661, 1666]; % line perpendicular to channel
%perpline = normdist;
%perpline = [perpline(2,:) - perpline(1,:)];
%perpline(2) = -perpline(2); % correcting for yaxis flipped in image coord
%[theta rho] = cart2pol(perpline(1), perpline(2));
mapscale = sum(round(sum(concentrationMap(1,:,:),3)));
%scalefactor = rho/mapscale;
%cart2pol 0 angle is at (1,0), and moves CCW, which matches degrees
% Hard coded channel width
channelWidth = 540.0316; % from new method
% previous hard code estimate is 1.219776924423975e+03
scalefactor = channelWidth/mapscale;

% Set rotation and center variables
%imangle = theta*180/pi; % degrees, CCW, this is coded to be 39.6606 deg CCW
imangle = 121.5889-90; % new one?
imscale = round(scalefactor); % what if don't round it?
imscale = scalefactor;
imcenter = mapcentroid;

% EXECUTE CONCENTRATION MAP
	concentrationMapRSC = RotateScaleCrop(concentrationMap, imscale, imangle, imcenter);

% Hard coded image center

imagecenter = [round(2101.18890250795);round(1565.36784139776)]; % hardcoded
% making Y of center smaller, moves the map DOWN
falseimage = zeros(3066,4086);
croppedMap = cropMap2Image(concentrationMapRSC, imageraw, imagecenter);




% Import cell data
calceinData = readtable('20200822_Live_Secondary.csv');
histogram(calceinData.Intensity_MeanIntensity_CropLive)
% from histogram learned that 0.02 is probably a good cutoff (left edge)
thresh = 0.02;

drug1C = croppedMap(:,:,1);
% find x and y locations of drug concentrations on image
x = calceinData.Location_Center_X;
y = calceinData.Location_Center_Y;

x = round(x);
y = round(y);

x = x + 1020; % because Kendall probably cropped before using image
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

%{
celldata = csvread('CellProfilerData_2019-07-09.csv',1,0);
cell_intensity_annexin = celldata(:,6); % annexin intensity
cell_intensity_calcein = celldata(:,7); % calcein intensity
cell_location = round(celldata(:, 40:41)); % XY reference
cell_intensity_ratio = cell_intensity_annexin./cell_intensity_calcein;
cell_N = length(cell_intensity_annexin); % number of cells

	% Exclude points that have 0 concentration values
	% First assign all points a score
	cell_concentration = zeros(cell_N, 3); % N by 3 matrix for concentration
	for i = 1:cell_N
		cell_concentration(i,:) = croppedMap(cell_location(i,2), cell_location(i,1),:); % indexing Y, X
	end
	cell_validity = find(sum(cell_concentration,2)==0); % points that are UNVALID
    cell_valid = find(sum(cell_concentration,2)~=0); % these are actually valid ones

% assemble a new matrix that eliminates 0 concentration points
% Structure [Object#, X, Y, C1, C2, C3, annexin, calcein, ratio]
cell_all = zeros(cell_N, 9);
	cell_all(:,1) = celldata(:,2); % object number
	cell_all(:,2:3) = cell_location; % location XY
	cell_all(:,4:6) = cell_concentration; % concentrations
	cell_all(:,7) = cell_intensity_annexin; % annexin intensity
	cell_all(:,8) = cell_intensity_calcein; % calcein intensity
	cell_all(:,9) = cell_intensity_ratio; % annexin/calcein ratio
cell_all(cell_validity,:) = []; % delete unvalid points




% Matrix of intensities in 2D spatial
intensity_image = zeros(size(croppedMap,1), size(croppedMap,2));
for i = 1:size(cell_all,1)
	intensity_image(cell_all(i,2),cell_all(i,3)) = cell_all(i,9);
end


[xq yq] = meshgrid(1:100:size(croppedMap,1),1:100:size(croppedMap,2));
vq = griddata(cell_all(:,2),cell_all(:,3),cell_all(:,9),xq,yq);
concentration_mask = round(sum(croppedMap,3));


%}

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
%{
% Figure 5. Image of annexin ratio intensity
figure(5); clf
imagesc(croppedMap(:,:,1))
hold on
for i = 1:cell_N
	plot(cell_location(i,1), cell_location(i,2),'.g')
end

for i = 1:length(cell_valid)
    plot(cell_location(cell_valid(i),1), cell_location(cell_valid(i),2), '.r')
end

hold off
colorbar
%}