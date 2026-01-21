function chamberSys = splittingChamber(f)
%SPLITTINGCHAMBER
% (Algebraic) Model of a splitting chamber with the input parameter f as
% the fraction of the first outlet flow of the overall inlet flow.

    eq = {"0 = f1*u1 - y1";...
        "0 = (1-f1)*u1 - y2"};
    
    chamberSys = sym2dmss(eq, 0);

    old = ["f1"];
    new  = [f];
    chamberSys = chamberSys.replaceSymbolicParameters(old, new);

end

