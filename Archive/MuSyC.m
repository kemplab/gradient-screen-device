function Ed = MuSyC(drug,varargin)
% Outputs drug effect based on 2D Hill surface
% Synergism determined by alpha and beta parameters
% Based on MuSyC system from Meyer et al. Cell Sys. 2019
% Supports multiple inputs of drugs up to 3

switch class(drug)
	case 'DrugSingle' % 1 drug Hill equation
		if nargin == 2
			% Extract parameters into "a"
			a(1) = drug.Emax;
			a(2) = drug.Hill;
			a(3) = drug.E0;
			a(4) = drug.EC50;
			Ed = MuSyC1(a, varargin{1});
		else
			disp('Not enough drug concentration inputs')
		end
		
	case 'DrugCombo' % 2 drug Hill equation
		if nargin == 3

			% Calculate E3 from beta
			E0 = drug.Drug1.E0;
            E1 = drug.Drug1.Emax;
			E2 = drug.Drug2.Emax;
			E3 = min(E1,E2)-drug.Beta*(E0-min(E1,E2));

			% Extract parameters into "a"
			a(1) = drug.Drug1.Emax;
			a(2) = drug.Drug1.Hill;
			a(3) = drug.Drug1.E0;
			a(4) = drug.Drug1.EC50;
			a(5) = drug.Drug2.Emax;
			a(6) = drug.Drug2.Hill;
			a(7) = drug.Drug2.EC50;
			a(8) = drug.Alpha(1);
			a(9) = drug.Alpha(2);
			a(10) = E3;
			a(11) = 2; % arbitrary number, unsure how to use functionally, but in paper
			a(12) = 2;
			Ed = MuSyC2(a, varargin{1}, varargin{2});
		else
			disp('Not enough drug concentration inputs')
		end
	
	case 'DrugTri' % 3 drug Hill equation
		if nargin == 4
			% Extract parameters into "a"
			%a(1) = drug.Emax;
			%a(2) = drug.Hill;
			%a(3) = drug.E0;
			%a(4) = drug.EC50;
			Ed = MuSyC3(drug, varargin{1}, varargin{2}, varargin{3});
		else
			disp('Not enough drug concentration inputs')
		end
	otherwise disp('Inputs must be drugSingle, drugCombo, or drugTri classes')
end

end

