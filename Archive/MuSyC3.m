function Ed = MuSyC3(tri,d1,d2,d3)
% 3D Hill equation, drug properties and concentration

h1 = tri.Drug1.Hill;
h2 = tri.Drug2.Hill;
h3 = tri.Drug3.Hill;
E0 = tri.Drug1.E0;
E1 = tri.Drug1.Emax;
E2 = tri.Drug2.Emax;
E3 = tri.Drug3.Emax;
C1 = tri.Drug1.EC50;
C2 = tri.Drug2.EC50;
C3 = tri.Drug3.EC50;
alpha = [tri.Combo12.Alpha, tri.Combo13.Alpha, tri.Combo23.Alpha];
beta = [tri.Combo12.Beta, tri.Combo13.Beta, tri.Combo23.Beta];
gamma = tri.Gamma;
r1 = 2; % arbitrary number
r2 = 2;
r3 = 2;
rn1 = (C1^h1)*r1; % updated from C1*r1 
rn2 = (C2^h2)*r2;
rn3 = (C3^h3)*r3;

conc1 = d1;
conc2 = d2;
conc3 = d3;

Ed = zeros(size(d1,1),size(d1,2),size(d1,3));


for i = 1:size(conc1,1)
        for j = 1:size(conc1,2)
                for k = 1:size(conc1,3) % Three nested loops in case of meshgrid

d1 = conc1(i,j,k);
d2 = conc2(i,j,k);
d3 = conc3(i,j,k);

        E12 = -beta(1)*min(E1,E2)+min(E1,E2);
        E13 = -beta(2)*min(E1,E3)+min(E1,E3);
        E23 = -beta(3)*min(E2,E3)+min(E2,E3);
        E123 = min([E12,E13,E23]); % how to infer the efficacy of drug 123 based on 3 pairs?

        Ev = [E0, E1, E2, E3, E12, E13, E23, E123];
        
        % U
        Y11 = - r1*d1^h1 - r2*d2^h2 - r3*d3^h3;
        Y12 = rn1;
        Y13 = rn2;
        Y14 = rn3;
        Y15 = 0;
        Y16 = 0;
        Y17 = 0;
        Y18 = 0;
        
        % A1
        Y21 = r1*d1^h1;
        Y22 = - rn1 - r2*alpha(1)*d2^h2 - r3*alpha(3)*d3^h3;
        Y23 = 0;
        Y24 = 0;
        Y25 = rn2;
        Y26 = rn3;
        Y27 = 0;
        Y28 = 0;
        
        % A2
        Y31 = r2*d2^h2;
        Y32 = 0;
        Y33 = - r1*alpha(2)*d1^h1 - rn2 - r3*alpha(5)*d3^h3;
        Y34 = 0;
        Y35 = rn1;
        Y36 = 0;
        Y37 = rn3;
        Y38 = 0;
        
        % A3
        Y41 = r3*d3^h3;
        Y42 = 0;
        Y43 = 0;
        Y44 = - r1*alpha(4)*d1^h1 - r2*alpha(6)*d2^h2 - rn3;
        Y45 = 0;
        Y46 = rn1;
        Y47 = rn2;
        Y48 = 0;
        
        % A12
        Y51 = 0;
        Y52 = r2*alpha(1)*d2^h2;
        Y53 = r1*alpha(2)*d1^h1;
        Y54 = 0;
        Y55 = - rn1 - rn2 - r3*gamma(1)*d3^h3;
        Y56 = 0;
        Y57 = 0;
        Y58 = rn3;
        
        % A13
        Y61 = 0;
        Y62 = r3*alpha(3)*d3^h3;
        Y63 = 0;
        Y64 = r1*alpha(4)*d1^h1;
        Y65 = 0;
        Y66 = - rn1 - r2*gamma(2)*d2^h2 - rn3;
        Y67 = 1;
        Y68 = rn2;
        
        % A23
        Y71 = 0;
        Y72 = 0;
        Y73 = r3*alpha(5)*d3^h3;
        Y74 = r2*alpha(6)*d2^h2;
        Y75 = 0;
        Y76 = 0;
        Y77 = - r1*gamma(3)*d1^h1 - rn2 - rn3;
        Y78 = rn1;
        
        % A123
        Y81 = 1;
        Y82 = 1;
        Y83 = 1;
        Y84 = 1;
        Y85 = 1;
        Y86 = 1;
        Y87 = 1;
        Y88 = 1;
        
        Y1 = [Y11, Y12, Y13, Y14, Y15, Y16, Y17, Y18];
        Y2 = [Y21, Y22, Y23, Y24, Y25, Y26, Y27, Y28];
        Y3 = [Y31, Y32, Y33, Y34, Y35, Y36, Y37, Y38];
        Y4 = [Y41, Y42, Y43, Y44, Y45, Y46, Y47, Y48];
        Y5 = [Y51, Y52, Y53, Y54, Y55, Y56, Y57, Y58];
        Y6 = [Y61, Y62, Y63, Y64, Y65, Y66, Y67, Y68];
        Y7 = [Y71, Y72, Y73, Y74, Y75, Y76, Y77, Y78];
        Y8 = [Y81, Y82, Y83, Y84, Y85, Y86, Y87, Y88];
        
        Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];
        
        b = [0;0;0;0;0;0;0;1];
        
        disp(Ev)
        disp(Y)
        disp(b)
        Ed(i,j,k) = Ev/Y*b;

end
end
end



end
% end
