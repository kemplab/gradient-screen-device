% create a slice of data that is a ternary plot slice simulating Kendall's device for James' data
data = readtable('data-hts-invert.csv');
conc1 = sort(unique(data.c1));
conc2 = sort(unique(data.c2));
conc3 = sort(unique(data.c3));

% extract 2-drug interaction tables
data12 = data(find(data.c3 == 0), :);
data13 = data(find(data.c2 == 0), :);
data23 = data(find(data.c1 == 0), :);

% create dose table as a test
data12_test = table;
data12_test.conc1 = data12.c1;
data12_test.conc2 = data12.c2;
data12_test.E12_bliss = data12.ECL1;

doseMatrix = doseTable2Matrix(data12_test);


% set max concentration of triangle
c1_max = 800;
c2_max = 752;
c3_max = 5.6;

% establish corners of triangle
P1 = [c1_max, 0, 0];
P2 = [0, c2_max, 0];
P3 = [0, 0, c3_max];

normal = cross(P1-P2, P1-P3);
normVec = normal./(norm(normal));

syms x y z
P = [x, y, z];

realdot = @(u,v) u*transpose(v);
planefunction = realdot(P-P1, normal);
%planefunction = dot(normal, P-P1);
zplane = solve(planefunction, z);

%x = linspace(0, c1_max, 20);
%y = -(c2_max/c1_max).*x + c2_max;
%z = eval(zplane);

x = linspace(0, 800, 10);
y = linspace(0, 752, 10);
[X,Y] = meshgrid(x,y);

% get rid of triangel bigger
X = X(:);
Y = Y(:);
idx = find(Y > -(c2_max/c1_max).*X + c2_max);
X(idx) = [];
Y(idx) = [];

x = X;
y = Y;
z = eval(zplane);
plot3(X,Y,z,'.')

xlabel('x1')
ylabel('x2')
zlabel('x3')

xnorm = x./c1_max;
ynorm = y./c2_max;
znorm = z./c3_max;

effect_test = rand(length(x), 1);
figure
addpath ternary2
tersurf(xnorm,ynorm,znorm,effect_test)


% interpolate the effect to points
effect = interp3(data.c1, data.c2, data.c3, data.ECL1, x, y, z);
figure
addpath ternary2
tersurf(xnorm,ynorm,znorm,effect)

% interp3 only works on data in meshgrid format
%Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
