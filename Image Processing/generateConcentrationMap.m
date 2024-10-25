function concentrationMap = generateConcentrationMap(mapGrid)
% Generates concentration map in image matrix format that can be manipulated
% Generates from [X Y C1 C2 C3] data

X = mapGrid(:,1);
Y = mapGrid(:,2);
C1 = mapGrid(:,3);
C2 = mapGrid(:,4);
C3 = mapGrid(:,5);

% remap XY spatial locations to a matrix
Xnew = normalizeCoordinates(X);
Ynew = normalizeCoordinates(Y);

% Accumulate into matrix
initialMapC1 = accumarray([Ynew, Xnew], C1, [max(Ynew), max(Xnew)]);
initialMapC2 = accumarray([Ynew, Xnew], C2, [max(Ynew), max(Xnew)]);
initialMapC3 = accumarray([Ynew, Xnew], C3, [max(Ynew), max(Xnew)]);
concentrationMap = zeros(size(initialMapC1,1),size(initialMapC1,2),3); % initialize size
concentrationMap(:,:,1) = flip(initialMapC1,1); % flip to correct image Y coordinates reversed
concentrationMap(:,:,2) = flip(initialMapC2,1);
concentrationMap(:,:,3) = flip(initialMapC3,1);