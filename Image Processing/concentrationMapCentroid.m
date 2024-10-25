function centroid = concentrationMapCentroid(concentrationMap)
% Outputs the centroid in [X Y] of a concentration map

concentrationMapBW = sum(concentrationMap,3);
concentrationMapBW(concentrationMapBW ~= 0) = 1;
s = regionprops(concentrationMapBW,'centroid'); % do object detection by contrast
if size(s.Centroid,1) > 1
	error('Multiple objects detected. Make sure map is correct.')
else
	centroid = s.Centroid;
end

centroid = round(centroid); % round to make real pixels