clear all;
close all;
clc;

% const
constant_phi = 0.968;
velocity_final = 340;
current_final = 0.48;
resistance_a = 1.63;
resistance_f = 658;
torque_simulation_range = -7:0.001:7;
num_samples = 43;
scaling_factor = 3;
interval_sampling = 1600;
interval_pre = 2000;
torque_final_values = [0.58; 0.45; 0.7];
voltage_a_values = [150; 50; 250];
error_tolerance = 0.05;
error_tolerance_max = 0.005;

% data filt
current_a_dcm_filtered = zeros(3, num_samples);
voltage_a_dcm_filtered = zeros(3, num_samples);
current_a_acm_filtered = zeros(3, num_samples);
voltage_an_acm_filtered = zeros(3, num_samples);
pf_acm_filtered = zeros(3, num_samples);
speed_n_filtered = zeros(3, num_samples);

measurement_files = {
    'Measurement_Data/Measurement_2023_12_14_15_44.mat',
    'Measurement_Data/Measurement_2023_12_14_15_57.mat',
    'Measurement_Data/Measurement_2023_12_14_16_1.mat'
};

for file_index = 1:length(measurement_files)
    load(measurement_files{file_index});
    for k = 1:num_samples
        range = interval_sampling*k*scaling_factor-interval_pre-1 : interval_sampling*k*scaling_factor-1;
        current_a_dcm_filtered(file_index, k) = mean(Ia_DCM(range));
        voltage_a_dcm_filtered(file_index, k) = mean(Va_DCM(range));
        current_a_acm_filtered(file_index, k) = mean(Ia_rms_ACM(range)); 
        voltage_an_acm_filtered(file_index, k) = mean(Van_rms_ACM(range)); 
        pf_acm_filtered(file_index, k) = mean(PF_ACM(range)); 
        speed_n_filtered(file_index, k) = mean(n(range));
    end
end

% model calc
torque_final_model = torque_final_values;
torque_simulation = torque_simulation_range;

voltage_a_model = voltage_a_values;

torque_drive_model = torque_simulation+torque_final_model;
current_a_model = torque_drive_model / constant_phi;
voltage_e_model = voltage_a_model - current_a_model*resistance_a;
speed_omega_model = voltage_e_model/constant_phi;
power_s_model = torque_simulation .* speed_omega_model;
speed_n_model = speed_omega_model/(2 * pi);

efficiency_model = (torque_simulation > 0).*(torque_drive_model > 0).*power_s_model./((voltage_a_model.*current_a_model)+(current_final*velocity_final));
efficiency_model = efficiency_model + (torque_simulation < 0).*(torque_drive_model < 0).*(abs(voltage_a_model.*current_a_model))./(abs(power_s_model)+(current_final*velocity_final));

% measure calc
speed_omega_meas = speed_n_filtered * 2*pi;

torque_final_meas = 0.8 * ((0.001638 * speed_omega_meas) + 0.4942);

torque_drive_meas = constant_phi*current_a_dcm_filtered;
torque_simulation_meas = torque_drive_meas-torque_final_meas;
power_s_meas = torque_simulation_meas .* speed_omega_meas;

efficiency_meas = (torque_simulation_meas > 0).*(torque_drive_meas > 0).*power_s_meas./((voltage_a_dcm_filtered.*current_a_dcm_filtered)+(current_final.*velocity_final));
efficiency_meas = efficiency_meas + (torque_simulation_meas < 0).*(torque_drive_meas < 0).*(abs(voltage_a_dcm_filtered.*current_a_dcm_filtered))./(abs(power_s_meas)+(current_final.*velocity_final));

% err calc
error_shaft = zeros(3, 1);
error_shaft_max = zeros(3, 1);
error_eff = zeros(3, 1);
error_eff_max = zeros(3, 1);

for i = 1:3
    error_shaft(i) = median(calculate_error_shaft(speed_omega_model(i,:), torque_simulation, speed_omega_meas(i,:), torque_simulation_meas(i,:), error_tolerance));
    error_shaft_max(i) = max(calculate_error_shaft(speed_omega_model(i,:), torque_simulation, speed_omega_meas(i,:), torque_simulation_meas(i,:), error_tolerance_max));
    error_eff(i) = median(calculate_error_efficiency(torque_simulation, efficiency_model(i,:), torque_simulation_meas(i,:), efficiency_meas(i,:), error_tolerance));
    error_eff_max(i) = max(calculate_error_efficiency(torque_simulation, efficiency_model(i,:), torque_simulation_meas(i,:), efficiency_meas(i,:), error_tolerance));
end

% plot
fig1 = figure(Units="inches", Position=[0 0 3.5 4]);
hold on;
model_colors = lines(3);

for k = 1:3
    plot(speed_omega_model(k,:), torque_simulation, 'LineWidth', 1.5, 'Color', model_colors(k,:));
    scatter(speed_omega_meas(k,:), torque_simulation_meas(k,:), 20, model_colors(k,:), 'filled', 'MarkerEdgeColor', 'k');
end

lgd = legend({'$V_a = 50\,V$', '$V_a = 50\,V$', ...
        '$V_a = 150\,V$', '$V_a = 150\,V$', ...
        '$V_a = 250\,V$', '$V_a = 250\,V$'}, ...
        'Interpreter', 'latex', 'Location', 'best');
lgd.ItemTokenSize(1) = 8;
xlabel('Shaft speed $\omega_m$ [rad/s]', 'Interpreter', 'latex');
ylabel('Shaft torque $T_s$ [Nm]', 'Interpreter', 'latex');
title('Shaft Speed vs Shaft Torque', 'Interpreter', 'latex');
grid on;
box on;
hold off;

% eff
fig2 = figure(Units="inches", Position=[0 0 3.5 4]);
hold on;

efficiency_meas = 100 * efficiency_meas;
efficiency_model = 100 * efficiency_model;

for k = 1:3
    plot(torque_simulation, efficiency_model(k,:), 'LineWidth', 1.5, 'Color', model_colors(k,:));
    scatter(torque_simulation_meas(k,:), efficiency_meas(k,:), 20, model_colors(k,:), 'filled', 'MarkerEdgeColor', 'k');
end

legend({'Model $V_a = 50\,V$', 'Measured $V_a = 50\,V$', ...
        'Model $V_a = 150\,V$', 'Measured $V_a = 150\,V$', ...
        'Model $V_a = 250\,V$', 'Measured $V_a = 250\,V$'}, ...
        'Interpreter', 'latex', 'Location', 'best');
xlabel('Shaft torque $T_s$ [Nm]', 'Interpreter', 'latex');
ylabel('Efficiency $\eta$ [\%]', 'Interpreter', 'latex');
title('Efficiency vs Shaft Torque', 'Interpreter', 'latex');
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