classdef DrugCombo
   properties
      Drug1
      Drug2
      Alpha
      Beta
   end
   methods
      function obj = DrugCombo(a,b)
         if nargin == 2
            obj.Drug1 = a;
            obj.Drug2 = b;
         end
      end
   end
end