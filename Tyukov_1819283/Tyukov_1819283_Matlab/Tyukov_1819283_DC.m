clear all;
close all;
clc;
% TODO: add constants from previous sections
% 5.5 kW Creusen 112L-4GM
k_phi = 0.968;
Vf = 340;
If = 0.48;
Ra = 1.63;
Rf = 658;
%%
% Calculate - Model
Tf_model = [0.58; 0.45; 0.7];
Ts_model = -7:0.001:7;

Va_model = [150; 50; 250];

Td_model = Ts_model+Tf_model;
Ia_model = Td_model / k_phi;
E_model = Va_model - Ia_model*Ra;
omega_m_model = E_model/k_phi;
Ps_model = Ts_model .* omega_m_model;
n_model = omega_m_model/(2 * pi);

eff_model = (Ts_model > 0).*(Td_model > 0).*Ps_model./((Va_model.*Ia_model)+(If*Vf));
eff_model = eff_model + (Ts_model < 0).*(Td_model < 0).*(abs(Va_model.*Ia_model))./(abs(Ps_model)+(If*Vf));

%%
% Calculate - Actual

Ns = 43;
s = 3; % time step

load('Measurement_Data/Measurement_2023_12_14_15_44.mat');
for k = 1:1:Ns
   Ia_DCM_filtered(1, k) = mean(Ia_DCM(1600*k*s-2000-1:1600*k*s-1));
   Va_DCM_filtered(1, k) = mean(Va_DCM(1600*k*s-2000-1:1600*k*s-1));
   Ia_ACM_filtered(1, k) = mean(Ia_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   Van_ACM_filtered(1, k) = mean(Van_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   PF_ACM_filtered(1, k) = mean(PF_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   n_filtered(1, k) = mean(n(1600*k*s-2000-1:1600*k*s-1));
end

load('Measurement_Data/Measurement_2023_12_14_15_57.mat');
for k = 1:1:Ns
   Ia_DCM_filtered(2, k) = mean(Ia_DCM(1600*k*s-2000-1:1600*k*s-1));
   Va_DCM_filtered(2, k) = mean(Va_DCM(1600*k*s-2000-1:1600*k*s-1));
   Ia_ACM_filtered(2, k) = mean(Ia_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   Van_ACM_filtered(2, k) = mean(Van_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   PF_ACM_filtered(2, k) = mean(PF_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   n_filtered(2, k) = mean(n(1600*k*s-2000-1:1600*k*s-1));
end

load('Measurement_Data/Measurement_2023_12_14_16_1.mat');
for k = 1:1:Ns
   Ia_DCM_filtered(3, k) = mean(Ia_DCM(1600*k*s-2000-1:1600*k*s-1));
   Va_DCM_filtered(3, k) = mean(Va_DCM(1600*k*s-2000-1:1600*k*s-1));
   Ia_ACM_filtered(3, k) = mean(Ia_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   Van_ACM_filtered(3, k) = mean(Van_rms_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   PF_ACM_filtered(3, k) = mean(PF_ACM(1600*k*s-2000-1:1600*k*s-1)); 
   n_filtered(3, k) = mean(n(1600*k*s-2000-1:1600*k*s-1));
end

% calculate omega_m continuous
omega_m_meas = n_filtered * 2*pi;

% continuous tfw_DCM
Tf_meas = 0.8 * ((0.001638 * omega_m_meas) + 0.4942);

Td_meas = k_phi*Ia_DCM_filtered;
Ts_meas = Td_meas-Tf_meas;
Ps_meas = Ts_meas .* omega_m_meas;

% combined eff of motor/gen
eff_meas = (Ts_meas > 0).*(Td_meas > 0).*Ps_meas./((Va_DCM_filtered.*Ia_DCM_filtered)+(If.*Vf));
eff_meas = eff_meas + (Ts_meas < 0).*(Td_meas < 0).*(abs(Va_DCM_filtered.*Ia_DCM_filtered))./(abs(Ps_meas)+(If.*Vf));

% Errors
% Inputs should have 1 row only
% Shaft torque - Speed
error_shaft = zeros(size(Ts_meas, 1), 1);
for i = 1:3
    error_shaft(i) = median(error_calc(omega_m_model(i,:), Ts_model, omega_m_meas(i,:), Ts_meas(i,:), 0.05));
end

error_shaft_max = zeros(size(Ts_meas, 1), 1);
for i = 1:3
    error_shaft_max(i) = max(error_calc(omega_m_model(i,:), Ts_model, omega_m_meas(i,:), Ts_meas(i,:), 0.005));
end

% Efficiency - Shaft torque
error_eff = zeros(size(eff_meas, 1), 1);
for i = 1:3
    error_eff(i) = median(error_calc(Ts_model, eff_model(i,:), Ts_meas(i,:), eff_meas(i,:), 0.05));
end

error_eff_max = zeros(size(eff_meas, 1), 1);
for i = 1:3
    error_eff_max(i) = max(error_calc(Ts_model, eff_model(i,:), Ts_meas(i,:), eff_meas(i,:), 0.05));
end

%%
figure;
hold on;
% Speed vs Shaft torque
col = ['r', 'g', 'b'];
for k=1:3
    plot(Ts_model, omega_m_model(k,:), Color=col(k));
end

for k=1:3
    plot(Ts_meas(k,:), omega_m_meas(k,:), 'o', Color=col(k));
end

legend('$V_a = 50 [V]$ ', '$V_a = 150 [V]$', '$V_a = 250 [V]$', 'interpreter', 'latex');
xlabel('Shaft speed $\omega_m$ [rad/s]', 'interpreter', 'latex');
ylabel('Shaft torque $T_s$ [Nm]', 'interpreter', 'latex');
hold off;

% Efficiency vs Speed
figure;
hold on;
for k=1:3
    plot(Ts_model, eff_model(k,:));
end
for k=1:3
    plot(Ts_meas(k,:), eff_meas(k,:), 'o');
end

legend('$V_a = 50 [V]$ ', '$V_a = 150 [V]$', '$V_a = 250 [V]$', 'interpreter', 'latex');
xlabel('Shaft speed $\omega_m$ [rad/s]', 'interpreter', 'latex');
ylabel('Shaft torque $T_s$ [Nm]', 'interpreter', 'latex');
hold off;