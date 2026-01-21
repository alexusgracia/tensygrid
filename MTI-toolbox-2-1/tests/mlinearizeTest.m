% Test header
function tests = mlinearizeTest
    tests = functiontests(localfunctions);
end

function firstTest(testCase)
    
   % Simulink File
    model = 'cstr_modelTest';
    
    % Parameters are also in the model workspace.
    
    % Time
    ts = 1e-2; %[h]
    t = 0:ts:1.8; % [h]
    
    % In workspace
    alpha = 30.828;
    delta = 3.556e-4;
    Hab = 4.2;
    beta = 86.688;
    gamma = 0.1;
    Hbc = -11;
    cin = 5.1e3;
    Tin = 104.9;
    Had = -41.85;
    E1 = 9758.3;
    k10 = 1.287e12;
    E2 = 8560;
    k20 = 9.043e6;
    
    % Initial Conditions
    cA0 = 0;
    cB0 = 0;
    Tr0 = 105;
    Tj0 = 105;
    x0 = [cA0;cB0;Tr0;Tj0];
    
    % Input function
    u = zeros(length(t),2);
    u(t>0.5 & t<1.5,1) = 2.5;
    u(t>=1.5,1) = 0.2;
    u(t>=1.5,2) = -1000;



   
    % Bounds go like [ca, cb, Tr, Tj, u, dq, outputs]
    up_bnd = [1000 750 116 116 2.5 0 1000 750 114 114];
    low_bnd = [0 0 104 104 0 -1000 0 0 104 104];

    % Nonlinear simulation with external input
    ds = Simulink.SimulationData.Dataset;
    ds = setElement(ds,1,timeseries(u(:,1),t));
    ds = setElement(ds,2,timeseries(u(:,2),t));

    simIn = Simulink.SimulationInput(model);

    simIn = simIn.setVariable('t',t);
    simIn = simIn.setVariable('alpha', alpha);
    simIn = simIn.setVariable('delta', delta);
    simIn = simIn.setVariable('Hab', Hab);
    simIn = simIn.setVariable('beta', beta);
    simIn = simIn.setVariable('gamma', gamma);
    simIn = simIn.setVariable('Hbc', Hbc);
    simIn = simIn.setVariable('cin', cin);
    simIn = simIn.setVariable('Tin', Tin);
    simIn = simIn.setVariable('Had', Had);
    simIn = simIn.setVariable('E1', E1);
    simIn = simIn.setVariable('k10', k10);
    simIn = simIn.setVariable('E2', E2);
    simIn = simIn.setVariable('k20', k20);
    simIn = simIn.setVariable('cA0', cA0);
    simIn = simIn.setVariable('cB0', cB0);
    simIn = simIn.setVariable('Tr0', Tr0);
    simIn = simIn.setVariable('Tj0', Tj0);
    simIn = simIn.setVariable('ts', ts);

    simIn = setExternalInput(simIn,ds);

    out = sim(simIn);  
    
    % Loop over different grid densities
    max_level = 6;
    min_level = 2;
    RMSE = zeros(max_level-min_level+1,4);
    %tmlin = zeros(1,max_level-min_level+1);
    for k=min_level:max_level
        
    % Sparse Projection Alg.
    level = k;
    max_order = 2;
    MLStic=tic;
    % State transion matrix and output matrix
    [msys] = mss.mlinearize(model,low_bnd,up_bnd,level,max_order);
    MLStoc=toc(MLStic);
    fprintf('%.2f sec - Multilinearisation\n', MLStoc)
    %tmlin(:,k-1) = MLStoc;
    
    % Perfom the simulation
    [~,tout,xout] = msim(msys,u,t,x0);

    % Continuos Simulation has time step defined by ODE solver, so
    % interpolate to match the dimensions
    x1_interpolated = interp1(tout(:,1), xout(:,1), t, 'linear')';
    x2_interpolated = interp1(tout(:,1), xout(:,2), t, 'linear')';
    x3_interpolated = interp1(tout(:,1), xout(:,3), t, 'linear')';
    x4_interpolated = interp1(tout(:,1), xout(:,4), t, 'linear')';

    x_int = [x1_interpolated, x2_interpolated, x3_interpolated, x4_interpolated];
    
    % STD of Error in trajectories
    RMSE(k-1,:) = std(out.yout-x_int);
    end
    
    % Check if all trajectory have pre-comupted error.
    comp = 0; % Number of comparisons. 
    for i=1:size(RMSE,1)  
        if all(RMSE(i,:)<[16 9 0.5 0.5])
            comp = comp+1;
        end
    end
    
    % If all trajectory errors are less than pre-computed error, then comp = max_level - min_level
    verifyEqual(testCase,comp-1,max_level-min_level);
 
end