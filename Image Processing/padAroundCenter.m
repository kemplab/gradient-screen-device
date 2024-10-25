function paddedImage = padAroundCenter(originalImage, centerPoint)
% Center Point in format [X, Y]

centerPoint = round(centerPoint); % turn this into integers

% Calculate true center, and find difference between true center and new center
centerX = (size(originalImage,2)+1)/2;
centerY = (size(originalImage,1)+1)/2;
%centerX = (size(originalImage,2))/2;
%centerY = (size(originalImage,1))/2;

diffX = centerX - centerPoint(1);
diffY = centerY - centerPoint(2);


% pad X first, if function to determine if places padding before or after
if diffX < 0
	paddedImage = padarray(originalImage,[0, abs(diffX*2)],'post');
else
	paddedImage = padarray(originalImage,[0, abs(diffX*2)],'pre');
end

% now pad Y
if diffY < 0
	paddedImage = padarray(paddedImage,[abs(diffY*2), 0],'post');
else
	paddedImage = padarray(paddedImage,[abs(diffY*2), 0],'pre');
end


end