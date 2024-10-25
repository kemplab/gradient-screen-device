function imageMap = remap2grid(X,Y,r,g,b)
% Converts vector points in X Y space into an image matrix in RGB space

Xnew = normalizeCoordinates(X); % normalize these to whole numbers
Ynew = normalizeCoordinates(Y);

% Accumulate into matrix
initialMapR = accumarray([Ynew, Xnew], r, [max(Ynew), max(Xnew)]);
initialMapG = accumarray([Ynew, Xnew], g, [max(Ynew), max(Xnew)]);
initialMapB = accumarray([Ynew, Xnew], b, [max(Ynew), max(Xnew)]);
imageMap = zeros(size(initialMapR,1),size(initialMapR,2),3); % initialize size
imageMap(:,:,1) = flip(initialMapR,1); % flip to correct image Y coordinates reversed
imageMap(:,:,2) = flip(initialMapG,1);
imageMap(:,:,3) = flip(initialMapB,1);
end