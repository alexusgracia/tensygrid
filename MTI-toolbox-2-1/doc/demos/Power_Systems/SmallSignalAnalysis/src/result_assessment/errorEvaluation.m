% --- Error Metrics Calculation Section ---
% Specify the time interval for error calculation
t_start = 2.5;   % or set externally as needed
t_end   = 2.505; % or set externally as needed

% --- GFM ---
[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gfmPQSimulink(:,1), timeDMSIM, gfmPQDMSIM(:,1), t_start, t_end);
errorEvaluation.P.multilinear.rmse = rmse;
errorEvaluation.P.multilinear.mae = mae;
errorEvaluation.P.multilinear.maxae = maxae;
errorEvaluation.P.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gfmPQSimulink(:,1), descriptorSimout.tsim, descriptorGFMpq(:,1), t_start, t_end);
errorEvaluation.P.descriptor.rmse = rmse;
errorEvaluation.P.descriptor.mae = mae;
errorEvaluation.P.descriptor.maxae = maxae;
errorEvaluation.P.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gfmPQSimulink(:,2), timeDMSIM, gfmPQDMSIM(:,2), t_start, t_end);
errorEvaluation.Q.multilinear.rmse = rmse;
errorEvaluation.Q.multilinear.mae = mae;
errorEvaluation.Q.multilinear.maxae = maxae;
errorEvaluation.Q.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gfmPQSimulink(:,2), descriptorSimout.tsim, descriptorGFMpq(:,2), t_start, t_end);
errorEvaluation.Q.descriptor.rmse = rmse;
errorEvaluation.Q.descriptor.mae = mae;
errorEvaluation.Q.descriptor.maxae = maxae;
errorEvaluation.Q.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, voltagePCCAmplitudeSimulink, timeDMSIM, voltageAmplitudeDMSIM, t_start, t_end);
errorEvaluation.V.multilinear.rmse = rmse;
errorEvaluation.V.multilinear.mae = mae;
errorEvaluation.V.multilinear.maxae = maxae;
errorEvaluation.V.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, voltagePCCAmplitudeSimulink, descriptorSimout.tsim, descriptorVload, t_start, t_end);
errorEvaluation.V.descriptor.rmse = rmse;
errorEvaluation.V.descriptor.mae = mae;
errorEvaluation.V.descriptor.maxae = maxae;
errorEvaluation.V.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gfmCurrentAmplitudeSimulink, timeDMSIM, gfmCurrentAmplitudeDMSIM, t_start, t_end);
errorEvaluation.I.multilinear.rmse = rmse;
errorEvaluation.I.multilinear.mae = mae;
errorEvaluation.I.multilinear.maxae = maxae;
errorEvaluation.I.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gfmCurrentAmplitudeSimulink, descriptorSimout.tsim, descriptorGFMcurrentAmplitude, t_start, t_end);
errorEvaluation.I.descriptor.rmse = rmse;
errorEvaluation.I.descriptor.mae = mae;
errorEvaluation.I.descriptor.maxae = maxae;
errorEvaluation.I.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, frequencySimulink, timeDMSIM, gfmFrequencyDMSIM, t_start, t_end);
errorEvaluation.f.multilinear.rmse = rmse;
errorEvaluation.f.multilinear.mae = mae;
errorEvaluation.f.multilinear.maxae = maxae;
errorEvaluation.f.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, frequencySimulink, descriptorSimout.tsim, descriptorGFMfrequency, t_start, t_end);
errorEvaluation.f.descriptor.rmse = rmse;
errorEvaluation.f.descriptor.mae = mae;
errorEvaluation.f.descriptor.maxae = maxae;
errorEvaluation.f.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, deltaAngleSimulink, timeDMSIM, deltaAngleDMSIM, t_start, t_end);
errorEvaluation.delta.multilinear.rmse = rmse;
errorEvaluation.delta.multilinear.mae = mae;
errorEvaluation.delta.multilinear.maxae = maxae;
errorEvaluation.delta.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, deltaAngleSimulink, descriptorSimout.tsim, descriptorGFMangle, t_start, t_end);
errorEvaluation.delta.descriptor.rmse = rmse;
errorEvaluation.delta.descriptor.mae = mae;
errorEvaluation.delta.descriptor.maxae = maxae;
errorEvaluation.delta.descriptor.nrmse = nrmse;

% --- GFL ---
[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gflPQSimulink(:,1), timeDMSIM, gflPQDMSIM(:,1), t_start, t_end);
errorEvaluation.P_gfl.multilinear.rmse = rmse;
errorEvaluation.P_gfl.multilinear.mae = mae;
errorEvaluation.P_gfl.multilinear.maxae = maxae;
errorEvaluation.P_gfl.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gflPQSimulink(:,1), descriptorSimout.tsim, descriptorGFLpq(:,1), t_start, t_end);
errorEvaluation.P_gfl.descriptor.rmse = rmse;
errorEvaluation.P_gfl.descriptor.mae = mae;
errorEvaluation.P_gfl.descriptor.maxae = maxae;
errorEvaluation.P_gfl.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gflPQSimulink(:,2), timeDMSIM, gflPQDMSIM(:,2), t_start, t_end);
errorEvaluation.Q_gfl.multilinear.rmse = rmse;
errorEvaluation.Q_gfl.multilinear.mae = mae;
errorEvaluation.Q_gfl.multilinear.maxae = maxae;
errorEvaluation.Q_gfl.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gflPQSimulink(:,2), descriptorSimout.tsim, descriptorGFLpq(:,2), t_start, t_end);
errorEvaluation.Q_gfl.descriptor.rmse = rmse;
errorEvaluation.Q_gfl.descriptor.mae = mae;
errorEvaluation.Q_gfl.descriptor.maxae = maxae;
errorEvaluation.Q_gfl.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, voltagePCCAmplitudeSimulink, timeDMSIM, voltageAmplitudeDMSIM, t_start, t_end);
errorEvaluation.V_gfl.multilinear.rmse = rmse;
errorEvaluation.V_gfl.multilinear.mae = mae;
errorEvaluation.V_gfl.multilinear.maxae = maxae;
errorEvaluation.V_gfl.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, voltagePCCAmplitudeSimulink, descriptorSimout.tsim, descriptorVload, t_start, t_end);
errorEvaluation.V_gfl.descriptor.rmse = rmse;
errorEvaluation.V_gfl.descriptor.mae = mae;
errorEvaluation.V_gfl.descriptor.maxae = maxae;
errorEvaluation.V_gfl.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gflCurrentAmplitudeSimulink, timeDMSIM, gflCurrentAmplitudeDMSIM, t_start, t_end);
errorEvaluation.I_gfl.multilinear.rmse = rmse;
errorEvaluation.I_gfl.multilinear.mae = mae;
errorEvaluation.I_gfl.multilinear.maxae = maxae;
errorEvaluation.I_gfl.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gflCurrentAmplitudeSimulink, descriptorSimout.tsim, descriptorGFLcurrentAmplitude, t_start, t_end);
errorEvaluation.I_gfl.descriptor.rmse = rmse;
errorEvaluation.I_gfl.descriptor.mae = mae;
errorEvaluation.I_gfl.descriptor.maxae = maxae;
errorEvaluation.I_gfl.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gflFrequencySimulink, timeDMSIM, gflFrequencyDMSIM, t_start, t_end);
errorEvaluation.f_gfl.multilinear.rmse = rmse;
errorEvaluation.f_gfl.multilinear.mae = mae;
errorEvaluation.f_gfl.multilinear.maxae = maxae;
errorEvaluation.f_gfl.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gflFrequencySimulink, descriptorSimout.tsim, descriptorGFLfrequency, t_start, t_end);
errorEvaluation.f_gfl.descriptor.rmse = rmse;
errorEvaluation.f_gfl.descriptor.mae = mae;
errorEvaluation.f_gfl.descriptor.maxae = maxae;
errorEvaluation.f_gfl.descriptor.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gfldeltaAngleSimulink, timeDMSIM, gfldeltaAngleDMSIM, t_start, t_end);
errorEvaluation.delta_gfl.multilinear.rmse = rmse;
errorEvaluation.delta_gfl.multilinear.mae = mae;
errorEvaluation.delta_gfl.multilinear.maxae = maxae;
errorEvaluation.delta_gfl.multilinear.nrmse = nrmse;

[rmse, mae, maxae, nrmse] = compute_errors(timeSimulink, gfldeltaAngleSimulink, descriptorSimout.tsim, descriptorGFLangle, t_start, t_end);
errorEvaluation.delta_gfl.descriptor.rmse = rmse;
errorEvaluation.delta_gfl.descriptor.mae = mae;
errorEvaluation.delta_gfl.descriptor.maxae = maxae;
errorEvaluation.delta_gfl.descriptor.nrmse = nrmse;

quantities = {'P', 'Q', 'V', 'I', 'f', 'delta'};
models = {'GFM Multilinear', 'GFM Descriptor'};

nrmse_mat = [
    errorEvaluation.P.multilinear.nrmse,  errorEvaluation.Q.multilinear.nrmse,  errorEvaluation.V.multilinear.nrmse,  errorEvaluation.I.multilinear.nrmse,  errorEvaluation.f.multilinear.nrmse,  errorEvaluation.delta.multilinear.nrmse;
    errorEvaluation.P.descriptor.nrmse,   errorEvaluation.Q.descriptor.nrmse,   errorEvaluation.V.descriptor.nrmse,   errorEvaluation.I.descriptor.nrmse,   errorEvaluation.f.descriptor.nrmse,   errorEvaluation.delta.descriptor.nrmse
];

maxae_mat = [
    errorEvaluation.P.multilinear.maxae,  errorEvaluation.Q.multilinear.maxae,  errorEvaluation.V.multilinear.maxae,  errorEvaluation.I.multilinear.maxae,  errorEvaluation.f.multilinear.maxae,  errorEvaluation.delta.multilinear.maxae;
    errorEvaluation.P.descriptor.maxae,   errorEvaluation.Q.descriptor.maxae,   errorEvaluation.V.descriptor.maxae,   errorEvaluation.I.descriptor.maxae,   errorEvaluation.f.descriptor.maxae,   errorEvaluation.delta.descriptor.maxae
];

T_nrmse = array2table(nrmse_mat, 'VariableNames', quantities, 'RowNames', models);
T_maxae = array2table(maxae_mat, 'VariableNames', quantities, 'RowNames', models);

disp('Normalized RMSE Table (GFM only):');
disp(T_nrmse);

disp('Maximum Absolute Error Table (GFM only):');
disp(T_maxae);


% Table for GFL converter: Maximum Absolute Error and Mean Absolute Error

quantities = {'P', 'Q', 'V', 'I', 'f', 'delta'};
models_gfl = {'GFL Multilinear', 'GFL Descriptor'};

mae_mat_gfl = [
    errorEvaluation.P_gfl.multilinear.mae,  errorEvaluation.Q_gfl.multilinear.mae,  errorEvaluation.V_gfl.multilinear.mae,  errorEvaluation.I_gfl.multilinear.mae,  errorEvaluation.f_gfl.multilinear.mae,  errorEvaluation.delta_gfl.multilinear.mae;
    errorEvaluation.P_gfl.descriptor.mae,   errorEvaluation.Q_gfl.descriptor.mae,   errorEvaluation.V_gfl.descriptor.mae,   errorEvaluation.I_gfl.descriptor.mae,   errorEvaluation.f_gfl.descriptor.mae,   errorEvaluation.delta_gfl.descriptor.mae
];

maxae_mat_gfl = [
    errorEvaluation.P_gfl.multilinear.maxae,  errorEvaluation.Q_gfl.multilinear.maxae,  errorEvaluation.V_gfl.multilinear.maxae,  errorEvaluation.I_gfl.multilinear.maxae,  errorEvaluation.f_gfl.multilinear.maxae,  errorEvaluation.delta_gfl.multilinear.maxae;
    errorEvaluation.P_gfl.descriptor.maxae,   errorEvaluation.Q_gfl.descriptor.maxae,   errorEvaluation.V_gfl.descriptor.maxae,   errorEvaluation.I_gfl.descriptor.maxae,   errorEvaluation.f_gfl.descriptor.maxae,   errorEvaluation.delta_gfl.descriptor.maxae
];

T_mae_gfl = array2table(mae_mat_gfl, 'VariableNames', quantities, 'RowNames', models_gfl);
T_maxae_gfl = array2table(maxae_mat_gfl, 'VariableNames', quantities, 'RowNames', models_gfl);

disp('Mean Absolute Error Table (GFL only):');
disp(T_mae_gfl);

disp('Maximum Absolute Error Table (GFL only):');
disp(T_maxae_gfl);


% Multiply P, Q, V, I columns by 100 (to get percent), f column by 1000 (to get mHz)
cols_percent = 1:4; % P, Q, V, I
col_freq = 5;       % f

% GFM tables
T_nrmse{:,cols_percent} = T_nrmse{:,cols_percent} * 100;
T_nrmse{:,col_freq} = T_nrmse{:,col_freq} * 1000;
T_maxae{:,cols_percent} = T_maxae{:,cols_percent} * 100;
T_maxae{:,col_freq} = T_maxae{:,col_freq} * 1000;

% GFL tables
T_mae_gfl{:,cols_percent} = T_mae_gfl{:,cols_percent} * 100;
T_mae_gfl{:,col_freq} = T_mae_gfl{:,col_freq} * 1000;
T_maxae_gfl{:,cols_percent} = T_maxae_gfl{:,cols_percent} * 100;
T_maxae_gfl{:,col_freq} = T_maxae_gfl{:,col_freq} * 1000;

% disp('Normalized RMSE Table (GFM only, % for P-Q-V-I, mHz for f):');
% disp(T_nrmse);
% 
% disp('Maximum Absolute Error Table (GFM only, % for P-Q-V-I, mHz for f):');
% disp(T_maxae);

disp('Maximum Absolute Error Table (GFL only, % for P-Q-V-I, mHz for f):');
disp(T_maxae_gfl);

disp('Mean Absolute Error Table (GFL only, % for P-Q-V-I, mHz for f):');
disp(T_mae_gfl);




% --- Error Metric Function ---
function [rmse, mae, maxae, nrmse] = compute_errors(reference_time, reference_signal, model_time, model_signal, t_start, t_end)
         % Ensure all incoming values are column vectors
    reference_time = reference_time(:);
    model_time = model_time(:);
    if isrow(reference_signal)
        reference_signal = reference_signal(:);
    end
    if isrow(model_signal)
        model_signal = model_signal(:);
    end   

    

    % Restrict to interval
    idx_ref = find(reference_time >= t_start & reference_time <= t_end);
    idx_mod = find(model_time >= t_start & model_time <= t_end);

    % % Ensure indices are within bounds
    % idx_ref = idx_ref(idx_ref <= size(reference_signal,1));
    % idx_mod = idx_mod(idx_mod <= size(model_signal,1));

    % Interpolate model to reference time base
    model_interp = interp1(model_time(idx_mod), model_signal(idx_mod), reference_time(idx_ref), 'linear', 'extrap');

    % Compute errors
    err = model_interp(:) - reference_signal(idx_ref);
    rmse = sqrt(mean(err.^2, 1));
    mae = mean(abs(err), 1);
    maxae = max(abs(err), [], 1);

    ref_range = max(reference_signal(idx_ref), [], 1) - min(reference_signal(idx_ref), [], 1);
    nrmse = (rmse ./ ref_range);
end