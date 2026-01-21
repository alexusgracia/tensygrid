function dYdt = cstr_system(t, Y, U)
    % Unpack the state variables
    T1 = Y(1);
    T2 = Y(3);
    T3 = Y(5);
    T4 = Y(7);
    Ca1 = Y(2);
    Ca2 = Y(4);
    Ca3 = Y(6);
    Ca4 = Y(8);
    
    % Inputs
    Qh1 = U(1);        % Heat input to the first reactor (KJ/h)
    Qh2 = U(2);        % Heat input to the second reactor (KJ/h)
    Qh3 = U(3);        % Heat input to the third reactor (KJ/h)
    Qh4 = U(4);        % Heat input to the fourth reactor (KJ/h)

    % Parameters
    F01 = 5;        % m3/h
    F1 = 35;        % m3/h
    F2 = 45;        % m3/h
    F3 = 33;        % m3/h
    F02 = 10;       % m3/h
    F03 = 8;        % m3/h
    F04 = 12;       % m3/h
    Fr1 = 20;       % m3/h
    Fr2 = 10;       % m3/h
    V1 = 1;         % m3/h
    V2 = 3;         % m3/h
    V3 = 4;         % m3/h
    V4 = 6;         % m3/h
    T01 = 300;      % K
    T02 = 300;      % K
    T03 = 300;      % K
    T04 = 300;      % K
    Ca01 = 4.0;     % kmol/m3
    Ca02 = 2.0;     % kmol/m3
    Ca03 = 3.0;     % kmol/m3
    Ca04 = 3.5;     % kmol/m3
    deltaH = [-5.0e4, -5.2e4, -5.0e4];  % KJ/kmol
    rho = 1000;     % kg/m3
    cp = 0.231;     % KJ/kg·K
    k0 = [3.0e6, 3.0e5, 3.0e5];  % h-1
    E = [5.0e4, 7.53e4, 7.53e4];  % KJ/kmol
    R = 8.314;      % KJ/kmol·K
 
    % Reaction rates
    r1 = k0.* exp(-E / (R * T1)) * Ca1;
    r2 = k0.* exp(-E / (R * T2)) * Ca2;
    r3 = k0.* exp(-E / (R * T3)) * Ca3;
    r4 = k0.* exp(-E / (R * T4)) * Ca4;

    % Temperature derivatives
    dT1dt = F01/V1 * (T01 - T1) + Fr1/V1 * (T2 - T1) + Fr2/V1 * (T4 - T1) - (sum(deltaH.*r1)/(rho * cp)) + Qh1 / (rho * cp * V1);
    dT2dt = F1/V2 * (T1 - T2) + F02/V2 * (T02 - T2) - (sum(deltaH.*r2)/(rho * cp)) + Qh2 / (rho * cp * V2);
    dT3dt = (F2 - Fr1)/V3 * (T2 - T3) + F03/V3 * (T03 - T3) - (sum(deltaH.*r3)/(rho * cp))+ Qh3 / (rho * cp * V3);
    dT4dt = F3/V4 * (T3 - T4) + F04/V4 * (T04 - T4) - (sum(deltaH.*r4)/(rho * cp)) + Qh4 / (rho * cp * V4);

    % Concentration derivatives
    dCa1dt = F01/V1 * (Ca01 - Ca1) + Fr1/V1 * (Ca2 - Ca1) + Fr2/V1 * (Ca4 - Ca1) - sum(r1);
    dCa2dt = F1/V2 * (Ca1 - Ca2) + F02/V2 * (Ca02 - Ca2) - sum(r2);
    dCa3dt = (F2 - Fr1)/V3 * (Ca2 - Ca3) + F03/V3 * (Ca03 - Ca3) - sum(r3);
    dCa4dt = F3/V4 * (Ca3 - Ca4) + F04/V4 * (Ca04 - Ca4) - sum(r4);

    % Combine derivatives into a column vector
    dYdt = [dT1dt; dCa1dt; dT2dt; dCa2dt; dT3dt; dCa3dt; dT4dt; dCa4dt];
end
