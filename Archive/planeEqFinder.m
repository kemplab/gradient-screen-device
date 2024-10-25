% create plane from box instead
x = linspace(0, 1, 10);
y = x;
z = x;

[X, Y, Z] = meshgrid(x,y,z);

X2 = X;
Y2 = Y;
Z2 = Z;

sumAll = X2 + Y2 + Z2;

idx1 = find(sumAll > 1+0.05 | sumAll < 1-0.05);


X2(idx1) = [];
Y2(idx1) = [];
Z2(idx1) = [];


plot3(X2(:),Y2(:),Z2(:), '.')