clear all;
close all;
clc;

% const
inductance_m = 0.3112;
inductance_s = 12.7e-3;
resistance_s = 1.84;
inductance_r = 6.3e-3;
resistance_r = 0.963;
current_ph = 5.6;
voltage_ph = 230;

% data load
num_samples = 43;
scaling_factor = 3;

measurement_files = {
    'Measurement_2024_1_11_15_10.mat',
    'Measurement_2024_1_11_15_15.mat',
    'Measurement_2024_1_11_15_19.mat'
};

current_a_dcm_filtered = zeros(3, num_samples);
voltage_a_dcm_filtered = zeros(3, num_samples);
current_a_acm_filtered = zeros(3, num_samples);
voltage_an_acm_filtered = zeros(3, num_samples);
pf_acm_filtered = zeros(3, num_samples);
speed_n_filtered = zeros(3, num_samples);

for file_idx = 1:length(measurement_files)
    load(measurement_files{file_idx});
    for k = 1:num_samples
        range = 1600*k*scaling_factor-2000-1 : 1600*k*scaling_factor-1;
        current_a_dcm_filtered(file_idx, k) = mean(Ia_DCM(range));
        voltage_a_dcm_filtered(file_idx, k) = mean(Va_DCM(range));
        current_a_acm_filtered(file_idx, k) = mean(Ia_rms_ACM(range)); 
        voltage_an_acm_filtered(file_idx, k) = mean(Van_rms_ACM(range)); 
        pf_acm_filtered(file_idx, k) = mean(PF_ACM(range)); 
        speed_n_filtered(file_idx, k) = mean(n(range));
    end
end

% const mod calc
constant_phi = 0.968;
speed_n_filtered(speed_n_filtered < -30) = NaN;
speed_omega_meas = speed_n_filtered * 2*pi;
torque_drive_dcm_meas = constant_phi*current_a_dcm_filtered;
torque_final_meas = 0.002401 * speed_omega_meas + 0.5246;
torque_drive_im_meas = torque_final_meas - torque_drive_dcm_meas;
torque_simulation_im_meas = 0.8*torque_final_meas - torque_drive_dcm_meas;

voltage_an_meas = voltage_an_acm_filtered;
current_a_meas = current_a_acm_filtered;
pf_meas = pf_acm_filtered;

efficiency_meas = (torque_simulation_im_meas > 0).*(torque_drive_im_meas > 0).*torque_simulation_im_meas.*speed_omega_meas./(3.*voltage_an_meas.*abs(current_a_meas).*pf_meas);
efficiency_meas = efficiency_meas + (torque_simulation_im_meas < 0).*(torque_drive_im_meas < 0).*(pf_meas < 0).*(3.*voltage_an_meas.*abs(current_a_meas).*abs(pf_meas))./(abs(torque_simulation_im_meas).*speed_omega_meas);

% mod calc
model_frequency = [40;25;10];
model_coefficient = polyfit(model_frequency, mean(voltage_an_acm_filtered(:,20:end), 2), 1);
coeff_b = model_coefficient(2);
coeff_a = model_coefficient(1);
voltage_an_model = coeff_a.*model_frequency + coeff_b;
slip_s = linspace(-1, 1, 1e4);
slip_s(slip_s==0) = NaN;

speed_omega_e_model = 2*pi.*model_frequency;
speed_omega_m_model = speed_omega_e_model.*(1-slip_s);
resistance_2 = resistance_r./slip_s;

% elec char
impedance_z1 = 1i.*speed_omega_e_model*inductance_s + resistance_s;
impedance_z2m = ((1i.*speed_omega_e_model*inductance_r + resistance_2).*(1i.*speed_omega_e_model*inductance_m))./(1i.*speed_omega_e_model*inductance_r + resistance_2 + 1i.*speed_omega_e_model*inductance_m);
impedance_ztot = impedance_z1 + impedance_z2m;
current_i2 = voltage_an_model./impedance_ztot.*(1i.*speed_omega_e_model*inductance_m)./(1i.*speed_omega_e_model*inductance_m + resistance_2 + 1i.*speed_omega_e_model*inductance_r);
current_a_model = voltage_an_model./impedance_ztot;

pf_model = cos(angle(current_a_model));

% mech char
torque_tf_tot_model = 0.002401 * speed_omega_m_model + 0.5246;
torque_tf_im_model = 0.2*torque_tf_tot_model;
torque_drive_im_model = 3*abs(current_i2).^2.*resistance_2./(speed_omega_e_model);
torque_simulation_im_model = torque_drive_im_model - torque_tf_im_model;

torque_simulation_im_model((torque_simulation_im_model < -7) | (torque_simulation_im_model > 7)) = NaN;
torque_simulation_im_meas((torque_simulation_im_meas < -7) | (torque_simulation_im_meas > 7)) = NaN;
speed_omega_m_model(isnan(torque_simulation_im_model)) = NaN;
speed_omega_meas(isnan(torque_simulation_im_meas)) = NaN;

% mot gen eff
efficiency_model = (torque_simulation_im_model > 0).*(torque_drive_im_model > 0).*(torque_simulation_im_model.*speed_omega_m_model)./(3.*voltage_an_model.*abs(current_a_model).*pf_model);
efficiency_model = efficiency_model + (torque_simulation_im_model < 0).*(torque_drive_im_model < 0).*(pf_model < 0).*(3.*voltage_an_model.*abs(current_a_model).*abs(pf_model))./(abs(torque_simulation_im_model).*abs(speed_omega_m_model));

% err coef
error_ts = zeros(3, 1);
error_ts_max = zeros(3, 1);
error_pf = zeros(3, 1);
error_pf_max = zeros(3, 1);
error_eff = zeros(3, 1);
error_eff_max = zeros(3, 1);

for i = 1:3
    error_ts(i) = median(calculate_error_shaft(speed_omega_m_model(i,:), torque_simulation_im_model(i,:), speed_omega_meas(i,:), torque_simulation_im_meas(i,:), 0.1));
    error_ts_max(i) = max(calculate_error_shaft(speed_omega_m_model(i,:), torque_simulation_im_model(i,:), speed_omega_meas(i,:), torque_simulation_im_meas(i,:), 0.06));
    error_pf(i) = median(calculate_error_efficiency(speed_omega_m_model(i,:), pf_model(i,:), speed_omega_meas(i,:), pf_meas(i,:), 0.04));
    error_pf_max(i) = max(calculate_error_efficiency(speed_omega_m_model(i,:), pf_model(i,:), speed_omega_meas(i,:), pf_meas(i,:), 0.01));
    error_eff(i) = median(calculate_error_efficiency(speed_omega_m_model(i,:), efficiency_model(i,:), speed_omega_meas(i,:), efficiency_meas(i,:), 0.005));
    error_eff_max(i) = max(calculate_error_efficiency(speed_omega_m_model(i,:), efficiency_model(i,:), speed_omega_meas(i,:), efficiency_meas(i,:), 0.002));
end

% plot
fig1 = figure('Units', 'inches', 'Position', [0 0 3.5 4]);
hold on;
color_palette = lines(3);

for k = 1:3
    scatter(speed_omega_meas(k,:), torque_simulation_im_meas(k,:), 20, color_palette(k,:), 'filled', 'MarkerEdgeColor', 'k');
    plot(speed_omega_m_model(k,:), torque_simulation_im_model(k,:), 'Color', color_palette(k,:), 'LineWidth', 1.5);
end

xlabel('Rotational Speed $\omega_m$ [rad/s]', 'Interpreter', 'latex');
ylabel('Torque on Shaft $T_s$ [Nm]', 'Interpreter', 'latex');
lgd = legend({'$f = 40$ Hz', '$f = 40$ Hz', ...
        '$f = 25$ Hz', '$f = 25$ Hz', ...
        '$f = 10$ Hz', '$f = 10$ Hz'}, ...
        'Interpreter', 'latex', 'Location', 'best');
lgd.ItemTokenSize(1) = 8;
grid on;
box on;
hold off;

% plot pf shs
fig2 = figure('Units', 'inches', 'Position', [0 0 3.5 4]);
hold on;

for k = 1:3
    scatter(speed_omega_meas(k,:), pf_meas(k,:), 20, color_palette(k,:), 'filled', 'MarkerEdgeColor', 'k');
    plot(speed_omega_m_model(k,:), pf_model(k,:), 'Color', color_palette(k,:), 'LineWidth', 1.5);
end

xlabel('Rotational Speed $\omega_m$ [rad/s]', 'Interpreter', 'latex');
ylabel('Power Factor', 'Interpreter', 'latex');
lgd = legend({'$f = 40$ Hz', '$f = 40$ Hz', ...
        '$f = 25$ Hz', '$f = 25$ Hz', ...
        '$f = 10$ Hz', '$f = 10$ Hz'}, ...
        'Interpreter', 'latex', 'Location', 'best');
lgd.ItemTokenSize(1) = 8;
grid on;
box on;
hold off;

% plot eff shs
efficiency_meas_percent = 100 * efficiency_meas;
efficiency_model_percent = 100 * efficiency_model;
fig3 = figure('Units', 'inches', 'Position', [0 0 3.5 4]);
hold on;

for k = 1:3
    scatter(speed_omega_meas(k,:), efficiency_meas_percent(k,:), 20, color_palette(k,:), 'filled', 'MarkerEdgeColor', 'k');
    plot(speed_omega_m_model(k,:), efficiency_model_percent(k,:), 'Color', color_palette(k,:), 'LineWidth', 1.5);
end

xlabel('Rotational Speed $\omega_m$ [rad/s]', 'Interpreter', 'latex');
ylabel('Efficiency $\eta$ [\%]', 'Interpreter', 'latex');
lgd = legend({'$f = 40$ Hz', '$f = 40$ Hz', ...
        '$f = 25$ Hz', '$f = 25$ Hz', ...
        '$f = 10$ Hz', '$f = 10$ Hz'}, ...
        'Interpreter', 'latex', 'Location', 'best');
lgd.ItemTokenSize(1) = 8;
grid on;
box on;
hold off;


% err sh
function error = calculate_error_shaft(simulated_speed, simulated_torque, measured_speed, measured_torque, tolerance)
    error = zeros(1, length(simulated_speed));
    error_idx = 1;

    measured_speed = measured_speed((measured_torque ~= 0) & (~isnan(measured_torque)));
    measured_torque = measured_torque((measured_torque ~= 0) & (~isnan(measured_torque)));

    for i = 1:length(measured_speed)
        x_idx = find(abs(simulated_speed - measured_speed(i)) <= tolerance);
        for j = 1:length(x_idx)
            err_value = 100*abs((simulated_torque(x_idx(j)) - measured_torque(i))/measured_torque(i));
            error(error_idx) = err_value;
            error_idx = error_idx + 1;
        end
    end
    error = error(1:error_idx-1);
end

% err eff
function error = calculate_error_efficiency(simulated_torque, simulated_eff, measured_torque, measured_eff, tolerance)
    error = zeros(1, length(simulated_torque));
    error_idx = 1;

    measured_torque = measured_torque((measured_eff ~= 0) & (~isnan(measured_eff)));
    measured_eff = measured_eff((measured_eff ~= 0) & (~isnan(measured_eff)));

    for i = 1:length(measured_torque)
        x_idx = find(abs(simulated_torque - measured_torque(i)) <= tolerance);
        for j = 1:length(x_idx)
            err_value = 100*abs((simulated_eff(x_idx(j)) - measured_eff(i))/measured_eff(i));
            error(error_idx) = err_value;
            error_idx = error_idx + 1;
        end
    end
    error = error(1:error_idx-1);
end