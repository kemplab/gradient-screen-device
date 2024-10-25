function distBetweenLines = lineNormDist(line1, line2)
% Finds the distance between two parallel lines
% Outputs [x1 y1; x2 y2] of the two points defining the line
% Lines are inputted in the format [X Y] for two points on the line
% If lines are not exactly parallel, the average of the slopes are used
% The intercepts are recalculated to fall in the midpoint of the old line

% line format
%[x1 y1; x2 y2]

% Generalizing the above, we have, d = |C1 − C2| / (A2 + B2)½
% For Ax + By + C = 0;
% Translating from y = B0 + B1X

coefficients1 = polyfit(line1(:,1)', line1(:,2)', 1); %1 deg polyfit
m1 = coefficients1(1);
b1 = coefficients1(2);

coefficients2 = polyfit(line2(:,1)', line2(:,2)', 1); %1 deg polyfit
m2 = coefficients2(1);
b2 = coefficients2(2);

m = (m1+m2)/2; % average the slopes

% use halfway point between points bounding line 1 to calculate norm line
startPoint = mean(line1,1); % still [X Y] format

% solve for the intercept of the other line
b3 = startPoint(2) + (1/m)*startPoint(1); % y + (1/m)x = b3

endPoint = zeros(1,2);
endPoint(1) = startPoint(1)+1; % random incremented X
endPoint(2) = polyval([-1/m, b3], endPoint(1)); % Y point from random X
line3 = [startPoint; endPoint];

% Solve equations for intersection of line 2 and 3
intersection = [m -1; -1/m -1]\[-b2; -b3]; % linear algebra AX = B
line4 = line3; % initialize line4 variable
line4(2,:) = intersection';


distBetweenLines = line4;

% Plot to show result
%{
figure
hold on
plot(line1(:,1), line1(:,2))
plot(line2(:,1), line2(:,2))
plot(line3(:,1), line3(:,2))
plot(line4(:,1), line4(:,2))
hold off
axis([0, 5, 0, 5])
%}