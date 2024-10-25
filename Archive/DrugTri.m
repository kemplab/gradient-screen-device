classdef DrugTri
   properties
      Drug1
      Drug2
      Drug3
      Combo12
      Combo13
      Combo23
      Gamma
   end
   methods
      function obj = DrugTri(a,b,c)
         if nargin == 3
            obj.Combo12 = a;
            obj.Combo13 = b;
            obj.Combo23 = c;
         end
      end
   end
end