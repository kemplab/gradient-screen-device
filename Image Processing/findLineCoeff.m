function [slope, intercept] = findLineCoeff(line)
% Finds the slope and intercept of a line from the drawline function
% line.Position

% line format
%[x1 y1; x2 y2]
x1 = line.Position(1,1);
y1 = line.Position(1,2);
x2 = line.Position(2,1);
y2 = line.Position(2,2);

m = (y2-y1)/(x2-x1);
% y = mx + b

b1 = y1 - m*x1;
b2 = y2 - m*x2;

if round(b1*1000)~=round(b2*1000) % 3 decimal precision
    error('Math is wrong.')
else
    b = b1;
end

slope = m;
intercept = b;



end