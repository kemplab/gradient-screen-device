function imageEdited = RotateScaleCrop(imageOriginal, imScale, rotAngle, centerPoint)
% Rotates, scales, and crops the original image per scaling factor and angle
% centerPoint needs to be an [x, y] coordinate

% Check first to make sure that point really is in the image
imHeight = size(imageOriginal,1);
imWidth = size(imageOriginal,2);
if centerPoint(1) <= imHeight && centerPoint(2) <= imWidth
else
	msg = 'Error: center point not within image bounds.'
	error(msg)
end

% PAD
imagePadded = padAroundCenter(imageOriginal, centerPoint);

% SCALE
imageScaled = imresize(imagePadded, imScale,'box');

% ROTATE
imageRotated = imrotate(imageScaled,rotAngle,'nearest','loose');

imageEdited = imageRotated;
end