% Data Generator function
% Generates simulated drug dose-response data based on inputs and model


%function doseResponseMatrix = dataGenerator(a, b)
% Inputs
%  - Single drug dose responses
% Outputs
%  - Expected viability


% Bliss data

% Two drug version first
% Initialize drugs
drugA = DrugSingle(0, 2, 1, 1);
drugB = DrugSingle(0, 2, 1, 1);

conc1 = linspace(0, 5, 100);
conc2 = linspace(0, 5, 100);

response1 = drugA.generateDoseResponse(conc1);
response2 = drugB.generateDoseResponse(conc2);


% Create dose table
doseFactorial = fullfact([length(conc1), length(conc2)]);
for i = 1:width(doseFactorial)
    doseList(:, i) = conc1(doseFactorial(:,i));
end
results = splitvars(table(doseList));
results = renamevars(results, ["doseList_1", "doseList_2"], ["conc1", "conc2"]);
results.E1 = zeros(height(results), 1);
results.E2 = zeros(height(results), 1);
results.E12_bliss = zeros(height(results), 1);
results.E1 = drugA.generateDoseResponse(results.conc1);
results.E2 = drugB.generateDoseResponse(results.conc2);
results.E12_bliss = results.E1 .* results.E2;

% visualize dose table 2D, collect into a grid
doseMatrix = doseTable2Matrix(results);

% select only diagonal samples, based on input doses
maxDose1 = 1;
maxDose2 = 2;
nRatios = 100; % number of ratiometric points along axis

% Axis direction is from max A to max B
diagDose1 = linspace(maxDose1, 0, nRatios);
diagDose2 = linspace(0, maxDose2, nRatios);
ratiometricDose = diagDose2 ./ (diagDose1 + diagDose2);
diagEffect = drugA.generateDoseResponse(diagDose1) .* drugB.generateDoseResponse(diagDose2);

figure
plot(ratiometricDose, diagEffect)
title('Bliss Null Diagonal Effect')
xlabel('Ratiometric Dose')
ylabel('Viability')
axis([0 1 0 1])

figure(1)
hold on
plot3(diagDose1, diagDose2, diagEffect, '.r')
hold off

%end % end of function
% -------------------