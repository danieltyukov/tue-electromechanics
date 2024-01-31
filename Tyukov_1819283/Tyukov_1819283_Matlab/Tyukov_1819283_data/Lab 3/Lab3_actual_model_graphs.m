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
td_IM_meas = 0.8*tf_tot_meas - td_DCM_meas;
ts_IM_meas = tf_tot_meas - td_DCM_meas;

Van_meas = Van_ACM_filtered;
Ia_meas = Ia_ACM_filtered;
PF_meas = PF_ACM_filtered;
eff_meas = (ts_IM_meas > 0).*(td_IM_meas > 0).*(ts_IM_meas.*omega_m_meas)./(3.*Van_meas.*abs(Ia_meas).*PF_meas);
% Generator
eff_meas = eff_meas + (ts_IM_meas < 0).*(td_IM_meas < 0).*(PF_meas < 0).*(3.*Van_meas.*abs(Ia_meas).*abs(PF_meas))./(abs(ts_IM_meas).*omega_m_meas);

%%
f_model = [40; 25; 10];
% b = 10, Van = 230 @f=50
% b = 10;
% a = (230-b)/50;

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


% Motor
eff_model = (ts_IM_model > 0).*(td_IM_model > 0).*(ts_IM_model.*omega_m_model)./(3.*Van_model.*abs(Ia_model).*PF_model);
% Generator
eff_model = eff_model + (ts_IM_model < 0).*(td_IM_model < 0).*(PF_model < 0).*(3.*Van_model.*abs(Ia_model).*abs(PF_model))./(abs(ts_IM_model).*abs(omega_m_model));
%%
%%plots
figure;
hold on;
col = ['r', 'g', 'b'];
for k=1:3
    meas(k) = plot(omega_m_meas(k,:), ts_IM_meas(k,:), 'o', Color=col(k));
    model(k) = plot(omega_m_model(k,:), ts_IM_model(k,:), Color=col(k));
end
yline(0);
title('Ts to omega');
hold off;

figure;
hold on;
for k=1:3
    meas(k) = plot(omega_m_meas(k,:), PF_meas(k,:), 'o', Color=col(k));
    model(k) = plot(omega_m_model(k,:), PF_model(k,:), Color=col(k));
end
title('Power factor');
hold off;

figure;
hold on;
for k=1:3
    meas(k) = plot(omega_m_meas(k,:), eff_meas(k,:), 'o', Color=col(k));
    model(k) = plot(omega_m_model(k,:), eff_model(k,:), Color=col(k));
end
title('Efficiency');
hold off;