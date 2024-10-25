classdef CellProps
	% Properties of cells detected in device, or groups of cells
   properties
      X % x spatial coordinates
      Y % y spatial coordinates
      C1 % concentration of drug 1
      C2 % concentration of drug 2
      C3 % concentration of drug 3
      Viability % average viability of cell calculated based on neighbors
      Live % live dead status of cell (boolean)
      Cluster % which cluster cell belongs in

   end
   methods
   end
end