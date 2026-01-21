% Load mss object 
clear all;
load("discretemodel/msysd.mat")

msys= mss(CPN1(msysd.F.U,msysd.F.phi),0.01);
%% Execution time evaluation for different n_concat values
% Note: Concat 2 doubles the model size, 5 quintuples etc..
% No new connections are made, so the model just gets bigger with repetetive information 
n_concat_values = [1 5 20 ];

% Number of repetition for averaging out times
num_repetitions = 3;

ts = 0.1; %  sample time
t = 0:ts:6; 
N = length(t);

% Heat inputs (base pattern)
Q1 = (1 + 2.*sin(0.2*pi*t)).*10^4;
Q2 = (2 + 4.*sin(0.2*pi*t)).*10^4;
Q3 = (2.5 + 5.*sin(0.2*pi*t)).*10^4;
Q4 = (1 + 2.*sin(0.2*pi*t)).*10^4;
U_base = [Q1;Q2;Q3;Q4];

% Base initial conditions
x0_base = [341.9213 3.3349 341.9161 3.0803 343.7129 3.1286 342.2733 3.3156]';

results_table = zeros(length(n_concat_values), 6); % n_concat, run1, run2, run3, mean, std
results_table_old = zeros(length(n_concat_values), 6); % n_concat, run1, run2, run3, mean, std

for i = 1:length(n_concat_values)
    n_concat = n_concat_values(i);
    execution_times = zeros(num_repetitions, 1);
    execution_times_old = zeros(num_repetitions, 1);
    
    % Create concatenated system once per n_concat value
    msyscc = conCatMss(msys, n_concat);
    x0 = repmat(x0_base, n_concat, 1);
    U = repmat(U_base, n_concat, 1);
    
    % Run simulation 3 times
    for rep = 1:num_repetitions
        msyscc.MEXaccelerate(true);
        % Equivalent call:
        % msyscc = MEXaccelerate(msyscc,true);
        tic_eval = tic;
        [ysim, tsim, xsim] = msim(msyscc, U', t', x0);
        execution_times(rep) = toc(tic_eval);
    end
    
    % Calculate statistics
    mean_time = mean(execution_times);
    std_time = std(execution_times);
    results_table(i, :) = [n_concat, execution_times', mean_time, std_time];

    % Change processXU to old version
    msyscc.MEXaccelerate(false);

    % Run simulation 3 times again
    for rep = 1:num_repetitions
        tic_eval = tic;
        [ysim, tsim, xsim] = msim(msyscc, U', t', x0);
        execution_times_old(rep) = toc(tic_eval);
    end
    
    % ReCalculate statistics
    mean_time_old = mean(execution_times_old);
    std_time_old = std(execution_times_old);
    

    % Store results
    results_table_old(i, :) = [n_concat, execution_times_old', mean_time_old, std_time_old];
    
end

summary = array2table(results_table,'VariableNames',["Num_Concat","1st run","2nd run","3rd run","Mean","Std"]);
summary_old = array2table(results_table_old,'VariableNames',["Num_Concat","1st run","2nd run","3rd run","Mean","Std"]);
disp("Old accumarray:")
disp(summary_old)
disp("New MTISIM:")
disp(summary)


function sysout = conCatMss(sys,n_concat)
    if(n_concat < 1)
        sysout = sys;
        return
    end
    sizephi = size(sys.F.phi);
    rank_old = sizephi(2);
    n = sys.n;
    m = sys.m;
    nm = n + m;
    FU = zeros(nm*(n_concat),rank_old*n_concat);
    Fphi = zeros(n*(n_concat) ,rank_old*n_concat);

    U = full(sys.F.U);
    phi = full(sys.F.phi);
    
    Fphi(1:n,1:rank_old) = phi;
    for i = 0:n_concat-1
        UStateStart = i*n+1;
        UStateEnd = ((i+1)*n);
        UInputStart = (n_concat)*n+1 + i*m;
        UInputEnd =  (n_concat)*n + i*m + m; 
        FU(UStateStart:UStateEnd,i*rank_old+1:(i+1)*rank_old) = U(1:n,:);
        FU(UInputStart:UInputEnd,i*rank_old+1:(i+1)*rank_old) = U(n+1:end,:);
        
        Fphi(i*n+1:(i+1)*n,i*rank_old+1:(i+1)*rank_old) = phi;
    end
    sysout = mss(CPN1(sparse(FU),sparse(Fphi)),sys.ts);
end