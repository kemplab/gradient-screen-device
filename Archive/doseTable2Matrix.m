% collect elements in a dose table into matrix form
% assume table has the following variables
%  - conc1
%  - conc2
%  - effect

% goal is to put it into a matrix, with conc1 and conc2 as the axes
function doseMatrix = doseTable2Matrix(doseTable)
conc1 = doseTable.conc1;
conc2 = doseTable.conc2;
effect = doseTable.E12_bliss;

% find unique doses
[doses1, indConc1, indDose1] = unique(conc1);
% doses = conc1(indConc), conc = doses(indDose);
[doses2, indConc2, indDose2] = unique(conc2);
[X, Y] = ndgrid(doses1, doses2); % can also use ndgrid, different orientation
Z = accumarray( [indDose1(:), indDose2(:)], effect, [], [], NaN ); % insert NaN if no value found

doseMatrix = Z;


surf(X,Y,Z)
xlabel('Drug 1')
ylabel('Drug 2')
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
view(0, 90) % top down view
colorbar