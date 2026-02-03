function  sysExchanger = wetHeatExchanger()   
%WETHEATEXCHANGER
% Creation of a hybrid Dmss Model of a cell of an Water-to-air cooler 
% including humidity saturation and dehumidification 
% with submodels of its water chamber, metal frame and air chamber.

    % Physical constants
    rhowater = 998;
    cw = 4200; %J/kgK
    
    rho = 1.204; %kg/m³
    cpd = 1006; %J/kgK
    cpv = 1860; %J/kgK
    h0 = 2501*10^3; %J/kg
    
    Pr = 0.7179;
    Sc = 0.6396;
    Le = Sc/Pr;
    
    Qcooler= -40000 / 4; %-10000 / 4; %W
    %Qcooler= -10000;
    ccooler = 477; %J/kgK 
    mcooler = 42.5 / 2; %42.5 / 4; %kg
    Rex = 3000;
    
    %curcuits = 14;
    curcuits = 28;
    %L = 2600 / 2;
    L = 1.301 / 2;
    %L = 1.301;
    di = 0.015; %m  (ca. 3/4 zoll)
    %di = 0.0129; %m (ca 3/8 zoll)
    Ai = di^2/4*pi();
    qv_water = 4290 / 3600 / 1000; %m³/s
    qm_water = qv_water * 998;
    mw = 998 *di^2/4 *pi()*L*curcuits*2;
    uwMax = qv_water /curcuits / (di^2/4*pi());
    alphaW = 780;
    DT = 0.5 ;%DT = 1;
    Aw = -Qcooler / (alphaW*DT);
    
    qv_air = 7000/3600 / 2; %m³/s
    %qv_air = 7000/3600;
    qm_air = rho* qv_air; % kg/s
    uaMax = 2.16; %m/s
    alphaA = 27;
    DT = 2;
    Aa = -Qcooler / (alphaA*DT);
    
    %% Water
    eqWater = ["mw*cw*xp1 = cw*u1*(u2-x1)+Aw*u4*(u3-x1)"];

    waterModel = sym2dmss(eqWater, 0);
    waterModel.stateName = ["outlet water temperature"];
    waterModel.inputName = ["water flow", "inlet water temperature", "plate temperature", "convective heattransfer coefficient water"];
    waterModel.H.phi.equality.c = subs(waterModel.H.phi.equality.c, "mw", mw);
    waterModel.H.phi.equality.c = subs(waterModel.H.phi.equality.c, "cw", cw);
    waterModel.H.phi.equality.c = double(subs(waterModel.H.phi.equality.c, "Aw", Aw));
    
    %% heat transfer water
    H = hyCPN1();
    
%     f = [0.4 0.8;...
%         0.2 0.4;...
%         0.75 0];
%     p =  [1 -1];

    f = [0.162 0.838;...
        0.398 0.571;...
        0.861, 0];
    p = [1 -1];
    
    H.F.input = [f(1:2,:)];
    H.F.algebraic = [f(3,:)];
    H.phi.equality = [p];
    convHeattransferWater = dmss(H,0);

    algebraicCurrentRange = [0 1];
    algebraicTargetRange = [101 828.5];
    inputCurrentRange = [0 1; 0 1];
    inputTargetRange = [[0.05 0.5]*curcuits*(di^2/4*pi())*998;...
                        0 50];
    convHeattransferWater = convHeattransferWater.normalize([],[], inputCurrentRange, inputTargetRange, algebraicCurrentRange, algebraicTargetRange);

    convHeattransferWater.algebraicName = ["convective heattransfer coefficient water"];
    convHeattransferWater.inputName = ["water flow", "inlet water temperature"];
    
    %% heat transfer air
    H = hyCPN1();
    
    f = [0.243 -0.2106;...
        -0.3096 0];
    p =  [0.1966 -0.1298];
    
    H.F.input = [f(1,:)];
    H.F.algebraic = [f(2,:)];
    H.phi.equality = p;
    
    convHeattransferAir = dmss(H,0);

    algebraicCurrentRange = [0 1];
    algebraicTargetRange = [26 26+33.2];
    %algebraicTargetRange = [13 26];
    inputCurrentRange = [0 1];
    inputTargetRange = [0.5 2.2]*qv_air/uaMax*rho;
    convHeattransferAir = convHeattransferAir.normalize([],[], inputCurrentRange, inputTargetRange, algebraicCurrentRange, algebraicTargetRange);

    convHeattransferAir.algebraicName = ["convective heattransfer coefficient air"];
    convHeattransferAir.inputName = ["air flow"];
    
    %% humidity saturation
    H = hyCPN1();
    H.F.input = [-0.016097 -0.015701 0 0;...
                0 0 0 1];
    H.F.algebraic = [0.503281 0 -1 0];
    H.F.boolean = [0 0 2 2];
    H.phi.equality = [-0.654361 0.325885 0 0];
    H.phi.inequality = [0 0 -1 -1];
    
    humiditySat = dmss(H,0);
    humiditySat.algebraicName = ["saturated humidity at plate temperature"];
    humiditySat.inputName = ["plate temperature" "inlet humidity"];
    humiditySat.booleanName = ["dehumidification activation"];
    
    %% Plate material
    eqPlate = ["0 = Aw*u1*(u2-y1) + Aa*(u4-y1)*u3 + u5*u3/(cpd*Le)*h0*Aa*(u6-u7)"];
    plateStorage = sym2dmss(eqPlate, 0);
    plateStorage.algebraicName = ["plate temperature"];
    plateStorage.inputName = ["convective heattransfer coefficient water", "outlet water temperature", "convective heattransfer coefficient air",...
        "outlet air temperature", "dehumidification activation",...
        "outlet humidity", "saturated humidity at plate temperature"];
    
    plateStorage = plateStorage.replaceSymbolicParameters(["Aw", "Aa", "Le", "cpd", "h0"],[Aw, Aa, Le, cpd, h0]);
    
    %% Air
    eqAir = ["0 = u1*(cpd*(u2-y1) + h0*(u3-y2) + cpv*(u2*u3-y1*y2)) + Aa*u4*(u5-y1) + u6*u4/(cpd*Le)*h0*Aa*(u7-y2)";...
        "0 = u1*(u3-y2) + u6*u4/(cpd*Le)*Aa*(u7-y2)"];
    
    airModel = sym2dmss(eqAir, 0);
    airModel.algebraicName = ["outlet air temperature", "outlet humidity"];
    airModel.inputName = ["air flow", "inlet air temperature", "inlet humidity", "convective heattransfer coefficient air",...
       "plate temperature", "dehumidification activation", "saturated humidity at plate temperature"];
    
    airModel = airModel.replaceSymbolicParameters(["cpd", "h0", "cpv", "Aa", "Le"],[cpd, h0, cpv, Aa, Le]);
    
    %% Connect everything
    sysExchanger = connect(plateStorage, convHeattransferAir, convHeattransferWater, waterModel, airModel, humiditySat);
    %clearvars -except sysExchanger
end