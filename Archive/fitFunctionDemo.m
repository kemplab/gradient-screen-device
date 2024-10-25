clear; clc; close all


data = importdata( 'OpioidHypnoticSynergy.txt' );
Propofol      = data.data(:,1);
Remifentanil  = data.data(:,2);
Algometry     = data.data(:,3);

Emax = 1;
Effect = @(IC50A, IC50B, alpha, n, x, y) ...
    Emax*( x/IC50A + y/IC50B + alpha*( x/IC50A )...
    .* ( y/IC50B ) ).^n ./(( x/IC50A + y/IC50B + ...
    alpha*( x/IC50A ) .* ( y/IC50B ) ).^n  + 1);

AlgometryEffect = fit( [Propofol, Remifentanil], Algometry, Effect, ...
    'StartPoint', [2, 10, 1, 0.8], ...
    'Lower', [-Inf, -Inf, -5, -Inf], ...
    'Robust', 'LAR' )

plot( AlgometryEffect, [Propofol, Remifentanil], Algometry )

