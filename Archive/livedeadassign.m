% analyzes calcein and PI levels using various methods to assign live-dead
%function live_1 = livedeadassign(data)

data = readtable('20210422_Dev12_concentrations.csv');

calcein = data.calceinMedian;
PI = data.PIMedian;

% multiple methods

% method 1: using calcein PI normalization equation
live_equation = calcein ./ (calcein + PI);
live_1 = logical(round(live_equation));

% method 2: using kmeans to separate into 2 clusters
clusters = kmeans([calcein, PI], 2);
live_2 = clusters==1;


%end