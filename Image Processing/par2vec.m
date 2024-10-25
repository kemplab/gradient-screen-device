function vec = par2vec(par)
% Transform the parametric form of 2D line (y=mx+b) to vector [x,y] form

% par is [m b]
% vec is [x y]

vec = zeros(1,2);
% start with x = 1;
m = par(1);
b = par(2);

x = 1;
y = m*x + b;

vec = [1 m];
vec = vec./norm(vec); % normalize to unit vector

end