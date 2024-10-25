clear;clc; close all


%% reducing the number of symbolic commands in the expression
d1 = linspace(0.001, 2, 20); % left side
d2 = linspace(2, 0.001, 20); % right side

y = real(test(d2, d1, -0.75));
plotaxis = linspace(0, 1, length(d1));
hold on
plot(plotaxis, y, '.-')
axis([0, 1, 0, 1])

y_additivity = (test(max(d1), 0, 0) - test(max(d2),0,0)).*plotaxis + test(max(d2),0,0);
hold on
plot(plotaxis, y_additivity, '.-')

% calculate the AUC from additivity
AUC_additivity = trapz(plotaxis, y_additivity)
AUC_data = trapz(plotaxis, y)
AUC = AUC_additivity - AUC_data





% calculate RBI at each point on the line?
E1 = test(d1, zeros(size(d1,1), size(d1,2)), 0);
E2 = test(d2, zeros(size(d2,1), size(d2,2)), 0);

y_lower_limit = 0.5 .* E1 .* E2;
plot(plotaxis, y_lower_limit, '--')

y_upper_limit = (E1.*E2) + (0.5 .* abs(min(E1, E2) - E1.*E2));
plot(plotaxis, y_upper_limit, '--')

y_bliss = E1.*E2;
plot(plotaxis, y_bliss, '--')

for i = 1:numel(y)
    DA_xy(i) = RBI(E1(i), E2(i), y(i));
end



% Functions ===================================================================
% actual function space

function y = test(d_1, d_2, a)
% a is interaction index from Greco et al. 1995
for i = 1:length(d_1);
ED50_1 = 1;
ED50_2 = 1;
E0 = 1;
Emax = 0;
h1 = 5;
h2 = 5;
d1 = d_1(i);
d2 = d_2(i);

% greco 1995 model, formulation of G L Drusano et al
syms E
termA = (d1 ./ (ED50_1 .* ((E - E0) ./ (Emax - E)).^ (1./h1)));
termB = (d2 ./ (ED50_2 .* ((E - E0) ./ (Emax - E)).^ (1./h2)));
termC = a .* d1 .* d2;
termD = ED50_1 .* ED50_2 .* ((E - E0) ./ (Emax - E)) .^ ((1 ./ (2 .* h1)) + (1 ./ (2 .* h2)));
equation = termA + termB + (termC ./ termD); % -1 to solve for equal to 0

test = vpasolve(equation == 1, E);
y(i) = eval(subs(test));

end % end of i loop
end % end of function
%test4 = vpasolve(equation == 1, E);
%y2 = eval(subs(test4));


% notes
% currently works for the vpasolve implementation of just d1
% trying to do this with the implementation for the entire equation



% Rescaled Bliss Independence
% Following method by Yeh's group
% Tekin et al. Ecology Letters. 2020
% Using a newly introduced framework to measure ecological stressor interactions

function [DA_xy, cat] = RBI(w_x, w_y, w_xy)

    %w_x = 0.70; % effect of stressor 1
    %w_y = 0.50; % effect of stressor 2
    w_xy_pred = w_x * w_y; % expected effect
    
    %w_xy = 0.6; % actual amount
    
    % Case of synergy
    if w_xy - w_xy_pred <= 0 % negative sign
        DA_xy = (w_xy - w_xy_pred) ./ abs(0 - w_xy_pred);
    
    elseif w_xy - w_xy_pred > 0
    
        % Case of antagonistic buffering
        if w_xy < min(w_x, w_y)
            DA_xy = (w_xy - w_xy_pred) ./ abs(min(w_x, w_y) - w_xy_pred);
              
        % Case of antagonistic suppression
        elseif w_xy > min(w_x, w_y)
            DA_xy = 1 + (w_xy - min(w_x, w_y)) ./ abs(1 - min(w_x, w_y));
        end
    
    end
    
    % Print result
    D_crit = 0.5; % threshold for criticality
    if abs(DA_xy) < D_crit
        cat = 'additivity';
    elseif DA_xy < 0
        cat = 'synergy';
    elseif DA_xy > 0 & DA_xy <= 1
        cat = 'antagonistic buffering';
    elseif DA_xy > 0 & DA_xy > 1
        cat = 'antagonistic suppression';
    else
        cat = 'Unknown'
    end
    
    end % end of function