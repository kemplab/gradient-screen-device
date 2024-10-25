classdef DrugSingle
   properties
      Name
      Hill
      EC50
      Emax
      E0
      R
   end
   methods
      % Class Constructor
      function obj = DrugSingle(a,b,c,d)
         if nargin == 4
            obj.Emax = a;
            obj.Hill = b;
            obj.E0 = c;
            obj.EC50 = d;         
         end
      end
      
      % Generate dose-response
      function doseResponse = generateDoseResponse(obj, conc)
          % conc is a vector of doses
          x = conc;
          Emax = obj.Emax;
          E0 = obj.E0;
          h = obj.Hill;
          C = obj.EC50;
          doseResponse = ((E0-Emax).*C.^h)./(C.^h + x.^h) + Emax;
      end
   end
end