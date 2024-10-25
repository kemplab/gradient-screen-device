clear;clc; close all

% rename variable names
data = readtable('xy_25um_c1c2c3.csv');
dev = renamevars(data, ["Var1", "Var2", "Var3", "Var4", "Var5"], ...
    ["x", "y", "c1", "c2", "c3"]);


% figure 1 concentration plot
figure
%[X, Y] = meshgrid(unique(dev.x), unique(dev.y));
%c1_grid = griddata(dev.x, dev.y, dev.c1, X, Y);
%c2_grid = griddata(dev.x, dev.y,dev.c2, X, Y);
%c3_grid = griddata(dev.x, dev.y,dev.c3, X, Y);
%pcolor(X,Y,c1_grid)
%shading interp
%colorbar
scatter(dev.x, dev.y, [], [dev.c1, dev.c2, dev.c3])


% figure 2 concentration 3d plot
figure
scatter3(dev.c1, dev.c2, dev.c3, 36, [dev.c1, dev.c2, dev.c3])
xlabel('Drug 1 concentration')
ylabel('Drug 2 concentration')
zlabel('Drug 3 concentration')


% ternary plot
figure
[ternX, ternY] = ternCoord(dev.c1, dev.c2, dev.c3);
scatter(ternX, ternY, [], [dev.c1, dev.c2, dev.c3])

% custom colormap matlab
figure
r = (0: 0.1 : 0.9)';
g = r.^1.8;
b = r.^2.1;
mymap = [r g b];
rgbplot(mymap)
hold on
colormap(mymap)
colorbar('Ticks', [])


% custom functions

% ternCoord: convert cartesian to ternary coordinate system
function [ternX, ternY] = ternCoord(a,b,c)
    ternX = 0.5.*(2.*b+c)./(a+b+c);
    ternY = (sqrt(3)./2).*c./(a+b+c);
end % end of function