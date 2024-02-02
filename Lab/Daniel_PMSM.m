clear all;
close all;
clc;

% const
resistance_r = 0.75;
inductance_l = 3.15e-3;
pole_pairs_p = 8;

% data calc
num_samples = 43;
time_step = 3;
current_a_dcm_filtered = zeros(3, num_samples);
voltage_a_dcm_filtered = zeros(3, num_samples);
current_a_acm_filtered = zeros(3, num_samples);
voltage_an_acm_filtered = zeros(3, num_samples);
pf_acm_filtered = zeros(3, num_samples);
speed_n_filtered = zeros(3, num_samples);

measurement_files = {
    'Measurement_Data/Measurement_2023_12_14_16_46.mat',
    'Measurement_Data/Measurement_2023_12_14_16_50.mat',
    'Measurement_Data/Measurement_2023_12_14_16_54.mat'
};

for file_idx = 1:length(measurement_files)
    load(measurement_files{file_idx});
    for k = 1:num_samples
        range = 1600*k*time_step-2000-1 : 1600*k*time_step-1;
        current_a_dcm_filtered(file_idx, k) = mean(Ia_DCM(range));
        voltage_a_dcm_filtered(file_idx, k) = mean(Va_DCM(range));
        current_a_acm_filtered(file_idx, k) = mean(Ia_rms_ACM(range)); 
        voltage_an_acm_filtered(file_idx, k) = mean(Van_rms_ACM(range)); 
        pf_acm_filtered(file_idx, k) = mean(PF_ACM(range)); 
        speed_n_filtered(file_idx, k) = mean(n(range));
    end
end

current_a_meas = current_a_acm_filtered;
speed_omega_meas = 2*pi*speed_n_filtered;
speed_omega_e_meas = speed_omega_meas*pole_pairs_p/2;

% model calc
beta_model = 0;
current_a_model = -5:0.005:5;
current_a_model = current_a_model(current_a_model ~= 0);

speed_n_model = [15; 25; 35];
speed_omega_model = 2*pi*speed_n_model;
speed_omega_e_model = speed_omega_model*pole_pairs_p/2;
flux_lambda_m = 0.133;

% elec char
voltage_ef_model = flux_lambda_m * speed_omega_e_model;
voltage_va_model = voltage_ef_model + current_a_model .* (resistance_r + 1i*speed_omega_e_model*inductance_l);

% tor
torque_tf_tot_model = 0.002401 * speed_omega_model + 0.5246;
torque_tf_pmsm_model = 0.2 * torque_tf_tot_model;
torque_td_model = 3*voltage_ef_model*current_a_model*cos(beta_model)./speed_omega_model;
torque_ts_model = torque_td_model - torque_tf_pmsm_model;

% pf eff
pf_model = cos(angle(voltage_va_model) - angle(current_a_model));
efficiency_model = (torque_ts_model > 0).*(torque_td_model > 0).*(torque_ts_model.*speed_omega_model)./abs(3*voltage_va_model.*current_a_model.*pf_model);
efficiency_model = efficiency_model + (torque_ts_model < 0).*(torque_td_model < 0).*abs(3*voltage_va_model.*current_a_model.*pf_model)./abs(torque_ts_model.*speed_omega_model);

% elec char
voltage_va_meas = voltage_an_acm_filtered.*exp(1i.*acos(pf_acm_filtered));
voltage_ef_meas = voltage_va_meas - current_a_meas.*(resistance_r + 1i.*speed_omega_e_meas*inductance_l);

% tor
beta_meas = angle(voltage_ef_meas);
torque_tf_tot_meas = 0.002401 * speed_omega_meas + 0.5246;
torque_tf_pmsm_meas = 0.2 * torque_tf_tot_meas;
torque_td_meas = 3*abs(voltage_ef_meas).*current_a_meas.*cos(beta_meas)./speed_omega_meas;
torque_ts_meas = torque_td_meas - torque_tf_pmsm_meas;

% pf eff
pf_meas = pf_acm_filtered;
pf_meas(abs(pf_meas) < 0.1) = NaN;
efficiency_meas = (torque_ts_meas > 0).*(torque_td_meas > 0).*(torque_ts_meas.*speed_omega_meas)./abs(3*voltage_va_meas.*current_a_meas.*pf_meas);
efficiency_meas = efficiency_meas + (torque_ts_meas < 0).*(torque_td_meas < 0).*abs(3*voltage_va_meas.*current_a_meas.*pf_meas)./abs(torque_ts_meas.*speed_omega_meas);

% err calc
error_pf = zeros(3, 1);
error_pf_max = zeros(3, 1);
error_eff = zeros(3, 1);
error_eff_max = zeros(3, 1);

for i = 1:3
    error_pf(i) = median(calculate_error_pf(torque_ts_model(i,:), pf_model(i,:), torque_ts_meas(i,:), pf_meas(i,:), 0.001));
    error_pf_max(i) = max(calculate_error_pf(torque_ts_model(i,:), pf_model(i,:), torque_ts_meas(i,:), pf_meas(i,:), 0.001));
    error_eff(i) = median(calculate_error_efficiency(torque_ts_model(i,:), efficiency_model(i,:), torque_ts_meas(i,:), efficiency_meas(i,:), 0.01));
    error_eff_max(i) = max(calculate_error_efficiency(torque_ts_model(i,:), efficiency_model(i,:), torque_ts_meas(i,:), efficiency_meas(i,:), 0.01));
end

% plot
fig1 = figure(Units="inches", Position=[0 0 3.5 4]);
hold on;
color_scheme = lines(3);

for k = 1:3
    scatter(torque_ts_meas(k,:), pf_meas(k,:), 20, color_scheme(k,:), 'filled', 'MarkerEdgeColor', 'k');
    plot(torque_ts_model(k,:), pf_model(k,:), 'Color', color_scheme(k,:), 'LineWidth', 1.5);
end

xlabel('Shaft torque $T_s$ [Nm]', 'Interpreter', 'latex');
ylabel('Power Factor', 'Interpreter', 'latex');
lgd = legend({'$\omega_m = 30\pi$ [rad/s]', '$\omega_m = 30\pi$ [rad/s]', ...
        '$\omega_m = 50\pi$ [rad/s]', '$\omega_m = 50\pi$ [rad/s]', ...
        '$\omega_m = 70\pi$ [rad/s]', '$\omega_m = 70\pi$ [rad/s]'}, ...
        'Interpreter', 'latex', 'Location', 'best');
lgd.ItemTokenSize(1) = 8;
grid on;
box on;
hold off;

fig2 = figure(Units="inches", Position=[0 0 3.5 4]);
hold on;
efficiency_meas_percent = 100 * efficiency_meas;
efficiency_model_percent = 100 * efficiency_model;

for k = 1:3
    scatter(torque_ts_meas(k,:), efficiency_meas_percent(k,:), 20, color_scheme(k,:), 'filled', 'MarkerEdgeColor', 'k');
    plot(torque_ts_model(k,:), efficiency_model_percent(k,:), 'Color', color_scheme(k,:), 'LineWidth', 1.5);
end

xlabel('Shaft torque $T_s$ [Nm]', 'Interpreter', 'latex');
ylabel('Efficiency $\eta$ [\%]', 'Interpreter', 'latex');
lgd = legend({'$\omega_m = 30\pi$ [rad/s]', '$\omega_m = 30\pi$ [rad/s]', ...
        '$\omega_m = 50\pi$ [rad/s]', '$\omega_m = 50\pi$ [rad/s]', ...
        '$\omega_m = 70\pi$ [rad/s]', '$\omega_m = 70\pi$ [rad/s]'}, ...
        'Interpreter', 'latex', 'Location', 'best');
lgd.ItemTokenSize(1) = 8;
grid on;
box on;
hold off;

% err calc pf
function error = calculate_error_pf(simulated_ts, simulated_pf, measured_ts, measured_pf, tolerance)
    error = zeros(1, length(simulated_ts));
    error_idx = 1;

    measured_ts = measured_ts((measured_pf ~= 0) & (~isnan(measured_pf)));
    measured_pf = measured_pf((measured_pf ~= 0) & (~isnan(measured_pf)));

    for i = 1:length(measured_ts)
        x_idx = find(abs(simulated_ts - measured_ts(i)) <= tolerance);
        for j = 1:length(x_idx)
            err_value = 100*abs((simulated_pf(x_idx(j)) - measured_pf(i))/measured_pf(i));
            error(error_idx) = err_value;
            error_idx = error_idx + 1;
        end
    end
    error = error(1:error_idx-1);
end

% err calc eff
function error = calculate_error_efficiency(simulated_ts, simulated_eff, measured_ts, measured_eff, tolerance)
    error = zeros(1, length(simulated_ts));
    error_idx = 1;

    measured_ts = measured_ts((measured_eff ~= 0) & (~isnan(measured_eff)));
    measured_eff = measured_eff((measured_eff ~= 0) & (~isnan(measured_eff)));

    for i = 1:length(measured_ts)
        x_idx = find(abs(simulated_ts - measured_ts(i)) <= tolerance);
        for j = 1:length(x_idx)
            err_value = 100*abs((simulated_eff(x_idx(j)) - measured_eff(i))/measured_eff(i));
            error(error_idx) = err_value;
            error_idx = error_idx + 1;
        end
    end
    error = error(1:error_idx-1);
end