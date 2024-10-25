function E = hill1D(a,x)
    Emax = a(1); % upper asymptote
    h = a(2); % hill coefficient
    E0 = a(3); % lower asymptote
    C = a(4); % EC50
    E = ((Emax-E0).*C.^h)./(C.^h + x.^h) + E0;
end