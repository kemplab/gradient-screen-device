function croppedMap = cropMap2Image(concentrationMap, immatrix, centroid)
% Takes a concentration map that has been properly scaled and rotated and centered
% Crops it to fit the exact size of a given image and the image center

W = size(immatrix,2); % width of image
H = size(immatrix,1); % height of image
x = centroid(1); % x location of centroid
y = centroid(2); % y location of centroid

% measurements of distance away from image center in 4 directions

%			difftop
%			   |
%   diffleft---O----diffright
%			   |
%			diffbot

difftop = -(1-y);
diffbot = H-y;
diffleft = -(1-x);
diffright = W-x;

% Get the centroid of the concentration map (to make calculations easier)
mapcentroid = concentrationMapCentroid(concentrationMap);
xmap = mapcentroid(1);
ymap = mapcentroid(2);

maptop = ymap - difftop;
mapbot = ymap + diffbot;
mapleft = xmap - diffleft;
mapright = xmap + diffright;

% verify that not going out of bounds TODO



croppedMap = concentrationMap(maptop:mapbot, mapleft:mapright, :);