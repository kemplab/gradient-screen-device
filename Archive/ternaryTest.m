clear;clc; clf
% Generate a grid of evenly spaced ternary data
a = linspace(0, 1, 20);
b = linspace(0, 1, 20);
[A, B] = meshgrid(a,b);
A2 = A(:);
B2 = B(:);
sumAB = A2+B2;
idx = find(sumAB>1);
A2(idx) = [];
B2(idx) = [];
C2 = ones(numel(A2),1) - A2 - B2;


% 2 degrees of freedom, c must depend on a+b
a = A2;
b = B2;
c = C2;
[ternX, ternY] = ternCoord(a,b,c);
plot(ternX, ternY, '.k')


% plot axes for reference
% Green - B=0
a_axis = linspace(0, 1, 100);
c_axis = 1-a_axis;
b_axis = zeros(1,numel(a_axis));
[ternX_B0, ternY_B0] = ternCoord(a_axis, b_axis, c_axis);
hold on
plot(ternX_B0, ternY_B0,'g-', 'LineWidth', 2)

% Blue - A=0
b_axis = linspace(0, 1, 100);
c_axis = 1-b_axis;
a_axis = zeros(1,numel(a_axis));
[ternX_A0, ternY_A0] = ternCoord(a_axis, b_axis, c_axis);
hold on
plot(ternX_A0, ternY_A0,'b-', 'LineWidth', 2)

% Yellow - C=0
b_axis = linspace(0, 1, 100);
a_axis = 1-b_axis;
c_axis = zeros(1,numel(a_axis));
[ternX_C0, ternY_C0] = ternCoord(a_axis, b_axis, c_axis);
hold on
plot(ternX_C0, ternY_C0,'y-', 'LineWidth', 2)

% Plot A max - blue
[ternX_Amax, ternY_Amax] = ternCoord(1, 0, 0);
hold on
plot(ternX_Amax, ternY_Amax, '.b', 'MarkerSize', 30)

% Plot B max - green
[ternX_Bmax, ternY_Bmax] = ternCoord(0, 1, 0);
hold on
plot(ternX_Bmax, ternY_Bmax, '.g', 'MarkerSize', 30)

% Plot C max - Yellow
[ternX_Cmax, ternY_Cmax] = ternCoord(0, 0, 1);
hold on
plot(ternX_Cmax, ternY_Cmax, '.y', 'MarkerSize', 30)

% Generate data for a fixed A/B ratio
hold on
a_ratio = linspace(0, 1, 100);
b_ratio = a_ratio.*2;
sumAB_ratio = a_ratio+b_ratio;
idx = find(sumAB_ratio>1);
a_ratio(idx) = [];
b_ratio(idx) = [];
c_ratio = 1-a_ratio-b_ratio;
[ternX_ratio, ternY_ratio] = ternCoord(a_ratio, b_ratio, c_ratio);
plot(ternX_ratio, ternY_ratio, '-r')
hold off

% Draw a line for 10% C
a_fixed = linspace(0, 1-0.1, 100);
c_fixed = 0.1.*ones(1,numel(a_fixed));
b_fixed = 1-c_fixed-a_fixed;
idx = find(a_fixed+b_fixed+c_fixed>1);
a_fixed(idx) = [];
b_fixed(idx) = [];
c_fixed(idx) = [];
[ternX_fixed, ternY_fixed] = ternCoord(a_fixed, b_fixed, c_fixed);
hold on
plot(ternX_fixed, ternY_fixed, 'y--','LineWidth',2)

function [ternX, ternY] = ternCoord(a,b,c)
ternX = 0.5.*(2.*b+c)./(a+b+c);
ternY = (sqrt(3)./2).*c./(a+b+c);
end
