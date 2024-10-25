% The four-parameter Hill equation
%
%                         d^h
% E = E0 + (Emax-E0) * ---------
%                      C^h + d^h

function [E0, Emax, C, h] = fit1D(x,E)
fun = @(a) hill1D(a,x) - E; % where a are the coefficients
x0 = [1, 1, 1, (max(x)-min(x))./2]; % initial guess
lb = [-inf, -inf, -inf, -inf]; % lower bound
ub = [inf, inf, inf, inf]; % upper bound
options = optimoptions('lsqnonlin','Display','None'); % disables text display
%options = optimoptions('lsqnonlin');
fitcoeffs = lsqnonlin(fun,x0,lb,ub,options); % solves
a = fitcoeffs;

Emax = a(1);
h = a(2);
E0 = a(3);
C = a(4);

figure
semilogx(x,E,'o')
x2 = logspace(log10(min(x)),log10(max(x)));
hold on
semilogx(x2,hill1D(a,x2),'-')
hold off
xlabel('dose')
ylabel('effect')
end