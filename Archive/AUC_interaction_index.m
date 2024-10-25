clear;clc;close all

% Objective
% Correlate the AUC to the interaction index alpha from Greco Model
% need to change this to a full-factorial simulation to get the effect
% instead of having nested for loops

% experimental design
%{
d1_max = 2;
d2_max = 2;
d1 = linspace(0.001, d1_max, 20); % left side
d2 = linspace(d2_max, 0.001, 20); % right side
h1 = [1, 2]; % h = 4 was too much for the system to handle apparently
h2 = [1, 2];
a = [-0.75, -0.25, 0, 0.25, 0.50, 0.75, 1, 2, 3];
%}


d1_max = 2;
d2_max = 2;
d1 = linspace(0.001, d1_max, 20); % left side
d2 = linspace(d2_max, 0.001, 20); % right side
h1 = [3]; % h = 4 was too much for the system to handle apparently
h2 = [1];
a = [0];
test_conditions = fullfact([numel(d1_max), numel(d2_max), numel(h1), numel(h2), numel(a)]);




for i = 1:height(test_conditions)
    y = greco(d1, d2, h1(test_conditions(i,3)), h2(test_conditions(i,4)), a(test_conditions(i,5)));
    AUC_ratio(i) = AUC_test(d1, d2, h1(test_conditions(i,3)), h2(test_conditions(i,4)), y);

    plotaxis = linspace(0, 1, length(d1));
    plot(plotaxis, y)
    hold on
end

axis([0, 1, 0, 1])
xlabel('ratiometric axis')
ylabel('effect')

figure
plot(a(test_conditions(:,5)), AUC_ratio,'.')
xlabel('interaction index')
ylabel('AUC')


% CUSTOM FUNCTIONS ============================================================

function y = greco(d_1, d_2, h1, h2, a)
    % a is interaction index from Greco et al. 1995
    y = zeros(1, length(d_1));
    for i = 1:length(d_1);
        ED50_1 = 1; % all input concentrations scaled to EC50
        ED50_2 = 1;
        E0 = 1;
        Emax = 0;
        %h1 = 5;
        %h2 = 5;
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
        y(i) = real(eval(subs(test))); % has a tendency to give imaginary numbers
    end % end of i loop

end % end of function


function AUC_ratio = AUC_test(d1, d2, h1, h2, y)
    % AUC stuff
    % determine line of additivity
    plotaxis = linspace(0, 1, length(d1));
    y_additivity = (greco(max(d1), 0, h1, h2, 0) - greco(max(d2),0, h1, h2, 0)).*plotaxis + greco(max(d2), 0, h1, h2, 0);

    % calculate the AUC from additivity
    AUC_additivity = trapz(plotaxis, y_additivity);
    AUC_data = trapz(plotaxis, y);
    AUC = AUC_additivity - AUC_data;
    AUC_ratio = AUC;
end % end of function