clear all;
close all;
clc;

% Parameters
% 3.0 kW Siemens 1AV3104A
Lm = 0.3112 % Magnetizing inductance
L1 = 12.7e-3 % Stator leakage inductance
R1 = 1.84 % Stator resistance
L2 = 6.3e-3 % Rotor leakage inductance
R2 = 0.963 % Rotor resistance
Iph = 5.6 % max current
Vph = 230 % max voltage

%%
% Calculate - Actual

Ns = 43;
s = 3; % time step

load('Measurement_Data/Measurement_2024_1_11_15_10.mat');
for k = 1:1:Ns
   Ia_DCM_filtered(1, k) = mean(Ia_DCM(1600*k*s-2000-1:1600*k*s-1));
   Va_DCM_filtered(1, k) = mean(Va_DCM(1600*k*s-2000-1:1600*k*s-1));
   Ia_ACM_filtered(1, k) = mean(Ia_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   Van_ACM_filtered(1, k) = mean(Van_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   PF_ACM_filtered(1, k) = mean(PF_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   n_filtered(1, k) = mean(n(1600*k*s-2000-1:1600*k*s-1));
end

load('Measurement_Data/Measurement_2024_1_11_15_15.mat');
for k = 1:1:Ns
   Ia_DCM_filtered(2, k) = mean(Ia_DCM(1600*k*s-2000-1:1600*k*s-1));
   Va_DCM_filtered(2, k) = mean(Va_DCM(1600*k*s-2000-1:1600*k*s-1));
   Ia_ACM_filtered(2, k) = mean(Ia_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   Van_ACM_filtered(2, k) = mean(Van_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   PF_ACM_filtered(2, k) = mean(PF_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   n_filtered(2, k) = mean(n(1600*k*s-2000-1:1600*k*s-1));
end

load('Measurement_Data/Measurement_2024_1_11_15_19.mat');
for k = 1:1:Ns
   Ia_DCM_filtered(3, k) = mean(Ia_DCM(1600*k*s-2000-1:1600*k*s-1));
   Va_DCM_filtered(3, k) = mean(Va_DCM(1600*k*s-2000-1:1600*k*s-1));
   Ia_ACM_filtered(3, k) = mean(Ia_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   Van_ACM_filtered(3, k) = mean(Van_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   PF_ACM_filtered(3, k) = mean(PF_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   n_filtered(3, k) = mean(n(1600*k*s-2000-1:1600*k*s-1));
end

% 5.5 kW Creusen 112L-4GM
k_phi = 0.968;
n_filtered(n_filtered < -30) = NaN;
% calculate omega_m continuous
omega_m_meas = n_filtered * 2*pi;
td_DCM_meas = k_phi*Ia_DCM_filtered;
tf_tot_meas = 0.002401 * omega_m_meas + 0.5246;
td_IM_meas = tf_tot_meas - td_DCM_meas;
ts_IM_meas = 0.8*tf_tot_meas - td_DCM_meas;

Van_meas = Van_ACM_filtered;
Ia_meas = Ia_ACM_filtered;
PF_meas = PF_ACM_filtered;

eff_meas = (ts_IM_meas > 0).*(td_IM_meas > 0).*(ts_IM_meas.*omega_m_meas)./(3.*Van_meas.*abs(Ia_meas).*PF_meas);
eff_meas = eff_meas + (ts_IM_meas < 0).*(td_IM_meas < 0).*(PF_meas < 0).*(3.*Van_meas.*abs(Ia_meas).*abs(PF_meas))./(abs(ts_IM_meas).*omega_m_meas);

%%
f_model = [40; 25; 10];
coeff = polyfit(f_model, mean(Van_ACM_filtered(:,20:end), 2), 1)
b = coeff(2);
a = coeff(1);
Van_model = a.*f_model + b;
s = linspace(-1, 1, 1e4);
s(s==0) = NaN;

omega_e_model = 2*pi.*f_model;
omega_m_model = omega_e_model.*(1-s);
R_2 = R2./s;
% electrical

Z1 = 1i.*omega_e_model*L1+R1;
Z2m = ((i.*omega_e_model*L2+R_2).*(i.*omega_e_model*Lm))./(i.*omega_e_model*L2+R_2+i.*omega_e_model*Lm);
Ztot = Z1+Z2m;
I2 = Van_model./Ztot.*(i.*omega_e_model*Lm)./(i.*omega_e_model*Lm + R_2 + i.*omega_e_model*L2);
Ia_model = Van_model./Ztot;

PF_model = cos(angle(Ia_model));

% mechanical

tf_tot_model = 0.002401 * omega_m_model + 0.5246;
tf_IM_model = 0.2*tf_tot_model;
td_IM_model = 3*abs(I2).^2.*R_2./(omega_e_model);
ts_IM_model = td_IM_model - tf_IM_model;

ts_IM_model((ts_IM_model < -7) | (ts_IM_model > 7)) = NaN;
ts_IM_meas((ts_IM_meas < -7) | (ts_IM_meas > 7)) = NaN;
omega_m_model(isnan(ts_IM_model)) = NaN;
omega_m_meas(isnan(ts_IM_meas)) = NaN;

% Motor
eff_model = (ts_IM_model > 0).*(td_IM_model > 0).*(ts_IM_model.*omega_m_model)./(3.*Van_model.*abs(Ia_model).*PF_model);
% Generator
eff_model = eff_model + (ts_IM_model < 0).*(td_IM_model < 0).*(PF_model < 0).*(3.*Van_model.*abs(Ia_model).*abs(PF_model))./(abs(ts_IM_model).*abs(omega_m_model));

% Shaft torque - Speed
error_ts = zeros(size(ts_IM_meas, 1), 1);
for i = 1:3
    error_ts(i) = median(error_calc(omega_m_model(i,:), ts_IM_model(i,:), omega_m_meas(i,:), ts_IM_meas(i,:), 0.1));
end
error_ts_max = zeros(size(ts_IM_meas, 1), 1);
for i = 1:3
    error_ts_max(i) = max(error_calc(omega_m_model(i,:), ts_IM_model(i,:), omega_m_meas(i,:), ts_IM_meas(i,:), 0.06));
end
% PF - Speed
error_pf = zeros(size(PF_meas, 1), 1);
for i = 1:3
    error_pf(i) = median(error_calc(omega_m_model(i,:), PF_model(i,:), omega_m_meas(i,:), PF_meas(i,:), 0.04));
end
error_pf_max = zeros(size(PF_meas, 1), 1);
for i = 1:3
    error_pf_max(i) = max(error_calc(omega_m_model(i,:), PF_model(i,:), omega_m_meas(i,:), PF_meas(i,:), 0.01));
end

% Efficiency - Speed
error_eff = zeros(size(eff_meas, 1), 1);
for i = 1:3
    error_eff(i) = median(error_calc(omega_m_model(i,:), eff_model(i,:), omega_m_meas(i,:), eff_meas(i,:), 0.005));
end
error_eff_max = zeros(size(eff_meas, 1), 1);
for i = 1:3
    error_eff_max(i) = max(error_calc(omega_m_model(i,:), eff_model(i,:), omega_m_meas(i,:), eff_meas(i,:), 0.002));
end
%%
%%plots
col = ['r', 'g', 'b'];
fig = figure(Units="inches");
fig.Position(3) = 3.5;
hold on;
for k=1:3
    plot(omega_m_meas(k,:), ts_IM_meas(k,:), 'o', Color=col(k));
end

for k=1:3
    plot(omega_m_model(k,:), ts_IM_model(k,:), Color=col(k));
end
legend('$f = 40$ [Hz]', '$f = 25$ [Hz]', '$f = 10$ [Hz]', 'interpreter', 'latex');
xlabel('Shaft speed $\omega_m$ [rad/s]', 'interpreter', 'latex');
ylabel('Shaft torque $T_s$ [Nm]', 'interpreter', 'latex');
hold off;

fig = figure(Units="inches");
fig.Position(3) = 3.5;
hold on;
for k=1:3
    plot(omega_m_meas(k,:), PF_meas(k,:), 'o', Color=col(k));
end

for k=1:3
    plot(omega_m_model(k,:), PF_model(k,:), Color=col(k));
end
legend('$f = 40$ [Hz]', '$f = 25$ [Hz]', '$f = 10$ [Hz]', 'interpreter', 'latex');
xlabel('Shaft torque $T_s$ [Nm]', 'Interpreter', 'latex');
ylabel('Power factor $\cos \beta$ [\%]', 'interpreter', 'latex');
hold off;

eff_meas = 100*eff_meas;
eff_model = 100*eff_model;

fig = figure(Units="inches");
fig.Position(3) = 3.5;
hold on;
for k=1:3
    plot(omega_m_meas(k,:), eff_meas(k,:), 'o', Color=col(k));
end

for k=1:3
    plot(omega_m_model(k,:), eff_model(k,:), Color=col(k));
end
legend('$f = 40$ [Hz]', '$f = 25$ [Hz]', '$f = 10$ [Hz]', 'interpreter', 'latex');
xlabel('Shaft torque $T_s$ [Nm]', 'Interpreter', 'latex');
ylabel('Efficiency $\eta$ [\%]', 'interpreter', 'latex');
hold off;

function error = error_calc(simulated_x, simulated_y, measured_x, measured_y, tolerance)
    % This function calculates the median absolute error between
    % simulations and measured data using a specified method for calculating errors.
    % Input: 
    %   simulated_x - x values of the simulated data
    %   simulated_y - y values of the simulated data
    %   measured_x - x values of the measured data
    %   measured_y - y values of the measured data
    %   tolerance - The tolerance within which x values are considered a match
    %   errorMethod - Function handle to the method used for calculating errors
    % Output: 
    %   median_absolute_error - The median absolute error between the datasets

    % Preallocate the errors array for maximum possible size
    error = zeros(1, size(simulated_x, 2));
    error_index = 1;

    measured_x = measured_x((measured_y ~= 0) & (~isnan(measured_y)));
    measured_y = measured_y((measured_y ~= 0) & (~isnan(measured_y)));
    % Loop through the measured array
    for i = 1:length(measured_x)
        % Find simulated_x values within the tolerance of the current measured_x value
        x_index = find(abs(simulated_x - measured_x(i)) <= tolerance);
        
        % Loop through x_index to calculate errors for all matching simulated_x values
        for j = 1:length(x_index)
            % Calculate error using the specified method
            err = 100*abs((simulated_y(x_index(j)) - measured_y(i))/measured_y(i));
            if (err > 100 | isnan(err))
                fprintf("%f at %f, %f at %f\n", simulated_y(x_index(j)), simulated_x(x_index(j)), measured_y(i), measured_x(i));
            end
            error(error_index) = err; % Store error
            error_index = error_index + 1;
        end
    end

    % Remove unused preallocated space
    error = error(1:error_index-1);
end
