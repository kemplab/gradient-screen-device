function MuSyC_coeffs = fit2DHill_update(y, x1, x2)

%{
MuSyC_coeffs definition
E1 = a(1); % Emax of drug 1
h1 = a(2); % Hill coefficient of drug 1
E0 = a(3); % E0 of drug 1 (and all)
C1 = a(4); % EC50 of drug 1

E2 = a(5); % Emax of drug 2
h2 = a(6); % Hill coefficient of drug 2
C2 = a(7); % EC50 of drug 2

alpha1 = a(8); % Assumes detailed balance
alpha2 = a(9);

E3 = a(10); % beta is used in the initialization of MuSyC to calculate E3
r1 = a(11);
r2 = a(12);
%}


modelComplexity = 2; % Only modelComplexity=2 is working right now

switch modelComplexity
	case 5 % rate of transition (r1, r2) >> 1
		% a1, a2, E3, E1, E2, C1, C2, h1, h2, E0, r1, r2
		x0 = [0.5, 1, 1, 0.1, 0.5, 1, 0.1, 1, 0, 0, 1, 1]; % initial guess
		lb = [0, 0, 0, 0, 0, 0, 0, -100, 0, 0, 0, 0]; % lower bound
		ub = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 2, 2]; % upper bound
		fun = @(a) MuSyC2(a,x1,x2) - y; % where a are the coefficients
	
	case 4 % System obeys detail balance
		% a2, E3, E1, E2, C1, C2, h1, h2, E0
		x0 = [0.5, 1, 1, 0.1, 0.5, 1, 0.1, 1, 0]; % initial guess
		lb = [0, 1, 0, 0, 0, 1, 0, -100, 0]; % lower bound
		ub = [100, 1, 100, 100, 100, 1, 100, 100, 100]; % upper bound
		fun = @(a) MuSyC2(a,x1,x2) - y; % where a are the coefficients
	
	case 3 % E0 is the minimally observed effect
		% a2, E3, E1, E2, C1, C2, h1, h2
		x0 = [0.5, 1, 1, 0.1, 0.5, 1, 0.1, 1, 0, 0, 2, 2]; % initial guess
		lb = [0, 1, 0, 0, 0, 1, 0, -100, 0, 0, 2, 2]; % lower bound
		ub = [100, 1, 100, 100, 100, 1, 100, 100, 100, 100, 2, 2]; % upper bound
		fun = @(a) MuSyC2(a,x1,x2) - y; % where a are the coefficients
	
	case 2 % h1, h2 are from single drug fits or 1 if single fits failed to converge
		% a2, E3, E1, E2, C1, C2
		
		% Find the 1D Hill fits for getting h1 and h2
		% Smallest dosing in C1
		C1smallest = min(unique(x1));
		C2smallest = min(unique(x2));

		indexd1const = find(x1==C1smallest); % this is used to fit response for drug2
		indexd2const = find(x2==C2smallest); % this is used to fit response for drug1

		drug1only = cellProps;
		drug1only.C1 = x1(indexd2const);
		drug1only.viability = y(indexd2const);
		drugFit1only = fit1DHill(drug1only);
		h1 = drugFit1only.hill;

		drug2only = cellProps;
		drug2only.C1 = x2(indexd1const);
		drug2only.viability = y(indexd1const);
		drugFit2only = fit1DHill(drug2only);
		h2 = drugFit2only.hill;

		E0 = max(y); % largest viability
		x0 = [0.1, min(y), drugFit1only.Emax, drugFit2only.Emax, drugFit1only.EC50, drugFit2only.EC50]; % initial guess
		lb = [0.000000001, 0, 0, 0, 0, 0]; % lower bound
		ub = [10000000, 100.1, 100.1, 100.1, max(x1), max(x2)]; % upper bound
		fun = @(a) MuSyC2(...
			[a(3), h1, E0, a(5), a(4), h2, a(6), a(1), a(1), a(2), 2, 2]...
		,x1,x2) - y; % where a are the coefficients

		%options = optimoptions('lsqnonlin','Display','None'); % disables text display
		options = optimoptions('lsqnonlin');
		fitcoeffs = lsqnonlin(fun,x0,lb,ub,options); % solves

		a = fitcoeffs;
		% for the order presented in line 79
        % a(1) = alpha2
        % a(2) = E3
        % a(3) = E1
        % a(4) = E2
        % a(5) = C1
        % a(6) = C2

		MuSyC_coeffs = [a(3), h1, E0, a(5), a(4), h2, a(6), a(1), a(1), a(2), 99, 99]; % r1 and r2 are arbitrary for now
		%{
		MuSyC_coeffs definition
		E1 = a(1); % Emax of drug 1
		h1 = a(2); % Hill coefficient of drug 1
		E0 = a(3); % E0 of drug 1 (and all)
		C1 = a(4); % EC50 of drug 1

		E2 = a(5); % Emax of drug 2
		h2 = a(6); % Hill coefficient of drug 2
		C2 = a(7); % EC50 of drug 2

		alpha1 = a(8); % Assumes detailed balance
		alpha2 = a(9);

		E3 = a(10); % beta is used in the initialization of MuSyC to calculate E3
		r1 = a(11);
		r2 = a(12);
		%}

	case 1 % C1, C2 are from single drug fits or the median concentration if single fits failed to converge
		% a2, E3, E1, E2
		x0 = [0.5, 1, 1, 0.1, 0.5, 1, 0.1, 1, 0, 0, 2, 2]; % initial guess
		lb = [0, 1, 0, 0, 0, 1, 0, -100, 0, 0, 2, 2]; % lower bound
		ub = [100, 1, 100, 100, 100, 1, 100, 100, 100, 100, 2, 2]; % upper bound
		fun = @(a) MuSyC2(a,x1,x2) - y; % where a are the coefficients
	
	case 0 % E1, E2 are assumed to be the maximally observed effect at maximum concentration of d1 and d2
		% a2, E3
		x0 = [0.5, 1, 1, 0.1, 0.5, 1, 0.1, 1, 0, 0, 2, 2]; % initial guess
		lb = [0, 1, 0, 0, 0, 1, 0, -100, 0, 0, 2, 2]; % lower bound
		ub = [100, 1, 100, 100, 100, 1, 100, 100, 100, 100, 2, 2]; % upper bound
		fun = @(a) MuSyC2(a,x1,x2) - y; % where a are the coefficients
end



end