function chamberSys = mixingChamber()
%MIXINGCHAMBER
% (Algebraic) Model of a mixing chamber for calculating mixed enthalpy,
% temperature and humidity.

    % Physical constants
    rhowater = 998;
    cw = 4200; %J/kgK
    
    rho = 1.204; %kg/m³
    cpd = 1006; %J/kgK
    cpv = 1860; %J/kgK
    h0 = 2501*10^3; %J/kg
    
    eq = {"0 = y1 - u1 - u4";...
        "0 = y1*y3 - u1*u3 - u4*u6";...
        "0 = y1*(cpd*y2+cpv*y2*y3+h0*y3) - u1*(cpd*u2+cpv*u2*u3+h0*u3) - u4*(cpd*u5+cpv*u5*u6+h0*u6)"};
    
    chamberSys = sym2dmss(eq, 0);

    old = ["cpd", "cpv", "h0"];
    new  = [cpd, cpv, h0];
    chamberSys = chamberSys.replaceSymbolicParameters(old, new);
end