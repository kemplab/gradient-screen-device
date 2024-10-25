% collect elements in a dose table into matrix form
% assume table has the following variables
%  - conc1
%  - conc2
%  - effect

% goal is to put it into a matrix, with conc1 and conc2 as the axes
function [c1_vec, c2_vec, doseMatrix] = drugResponseVec2Mat(c1, c2, effect)



% find unique doses
[doses1, indConc1, indDose1] = unique(c1);
% doses = conc1(indConc), conc = doses(indDose);
[doses2, indConc2, indDose2] = unique(c2);
[c1_vec, c2_vec] = ndgrid(doses1, doses2); % can also use ndgrid, different orientation
Z = accumarray( [indDose1(:), indDose2(:)], effect, [], [], NaN ); % insert NaN if no value found

doseMatrix = flipud(Z);
c1_vec = flipud(c1_vec);
c2_vec = flipud(c2_vec);

surf(c1_vec,c2_vec,doseMatrix)
xlabel('Drug 1')
ylabel('Drug 2')
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
view(0, 90) % top down view
colorbar