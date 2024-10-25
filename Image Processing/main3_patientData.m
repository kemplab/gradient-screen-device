filename = '20210422_Dev17';
Calcein = readtable([filename, 'Calcein', '.csv']);
PI = readtable([filename, 'PI', '.csv']);
x = Calcein{:,21};
y = Calcein{:,22};

plot(x,y,'.')

% After this, all copied from main script


lines = zeros(6,2); % 6 lines, slope and intercept
h = msgbox('Trace all the lines bordering channels using the cursor (n=6)');
waitfor(h)

for i = 1:6
    rawlines(i) = drawline('LineWidth',1,'Color','g','Label',num2str(i));
    [lines(i,1), lines(i,2)] = findLineCoeff(rawlines(i));
end
uiwait(msgbox('This message will pause execution until you click OK. Double check your lines'));

% index pairs of parallel lines
parallellines = zeros(3,2); % 3 pairs, 2 indices of lines that are parallel
parscore = zeros(6,6);
for i = 1:6 % all lines
    for j = 1:6 % all pairs give a parallel score
        if i==j % make sure line doesn't compare itself
            parscore(i,j) = NaN;
        else
            parscore(i,j) = (lines(i,1)-lines(j,1));
        end
    end
end
%parscore(parscore < 1) = NaN; % if it changes sign, definitely not a line
%parscore(parscore > 150) = NaN; % if it's too large, there's too much of a difference, changed from 100
% always a reciprocal pair, let's only take 1 (hence < 1 and not < 0)
%[pairs1, pairs2] = find(parscore<20); % changed this threshold from 1.5, maybe can't be a ratio
%if length(pairs1)~=3
%    error('Program did not find 3 pairs of parallel lines. Please try again.')
%end
%parallellines(:,1) = pairs1;
%parallellines(:,2) = pairs2;

parallellines = [1,2; 3,4; 5,6];




%% CENTERLINES AND FINDING CENTER

% Find centerlines
centerline = zeros(3,2); % 3 lines, 2 parametric points (m,b)
for i = 1:3 % three pairs of lines
    centerline(i,1) = mean(...
        [lines(parallellines(i,1),1), lines(parallellines(i,2),1)]); % m
    centerline(i,2) = mean(...
        [lines(parallellines(i,1),2), lines(parallellines(i,2),2)]); % b
end


% Plot centerlines
hold on
X = 1:max(x);
for i = 1:3
    m = centerline(i,1);
    b = centerline(i,2);
    Y = m.*X + b;
    plot(X,Y,'-r');
end
hold off



% Find line intersection
m = centerline(:,1);
b = centerline(:,2);
intersection = [m, -ones(3,1)]\-b; % linear algebra AX = B
intersection2 = [lines(:,1), -ones(6,1)]\-lines(:,2); % all the lines instead

hold on
plot(intersection(1), intersection(2), 'ro')
plot(intersection2(1), intersection2(2), 'mo')
hold off

% intersection 1 is 2993, 3014
% intersection 2 is 2995, 3013


% Set the lines 1 and 2 to be the "top" channel, set all rotations based on that
topChannel = 1;
% Color the corresponding lines the color blue
set(rawlines(parallellines(topChannel,1)),'Color','b')
set(rawlines(parallellines(topChannel,2)),'Color','b')



%% FIND THE WIDTH OF EACH CHANNEL IN PIXELS
% Uses lineNormDist function
for i = 1:3
    normline = lineNormDist(rawlines(parallellines(i,1)).Position,...
        rawlines(parallellines(i,2)).Position);
    normline = [normline(2,:) - normline(1,:)];
    [~, rho] = cart2pol(normline(1),normline(2));
    channelWidth(i) = rho;
end
meanChannelWidth = mean(channelWidth);
% this turns out to be 1.16300343057825e+03




%% FIND THE ROTATION OF THE IMAGE IN RELATION TO THE TOP CHANNEL
% Use the centerline
%vectorofcenterline = par2vec(centerline(topChannel,:));
%[theta, rho] = cart2pol(-vectorofcenterline(1),vectorofcenterline(2));
%angleRotate = theta*180/pi; % This turns out to be 128.6834 degrees

angleRotate = atand(-1/centerline(1)); % degrees to turn counterclockwise from y-axis


% need to rotate by angleRotate
% need to scale by channelWidth
% center by intersection1(x,y)
