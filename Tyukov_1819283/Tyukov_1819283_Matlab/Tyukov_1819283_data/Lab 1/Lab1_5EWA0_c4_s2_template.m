%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         
%   Script name:            Lab1_5EWA0_c4_s2_template.m                      
%   Comments:               Template script for lab 1 question 4.2
%   Last checked/changed:   October 2023 by M. Kleijer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

% TODO: save this script in your own work folder. Remove template from the
% name and add your own name.

%% Import data from LabVIEW
% Datafile should be in same folder as this script.

load('Measurement_Data/Measurement_2023_12_6_17_24.mat');

% TODO: add constants from previous sections
k_phi = 0.968;
Vf = 327.6;
If = 0.46;
Ra = 1.63;
Rf = 658;

%% Filter data 

% TODO: input the number of steps (Ns) and the step time (Ts) used in the
% torque-sequencer:
Ns = 43;
Ts = 3;


for i = 1:1:Ns
   Ia_DCM_filtered(i) = mean(Ia_DCM(1600*i*Ts-2000-1:1600*i*Ts-1));
   Va_DCM_filtered(i) = mean(Va_DCM(1600*i*Ts-2000-1:1600*i*Ts-1));
   Ia_ACM_filtered(i) = mean(Ia_rms_ACM(1600*i*Ts-2000-1:1600*i*Ts-1)); 
   Van_ACM_filtered(i) = mean(Van_rms_ACM(1600*i*Ts-2000-1:1600*i*Ts-1)); 
   PF_ACM_filtered(i) = mean(PF_ACM(1600*i*Ts-2000-1:1600*i*Ts-1)); 
   n_filtered(i) = mean(n(1600*i*Ts-2000-1:1600*i*Ts-1));
end
clear T_seq;                                % Not used in the 5EWA0 labs
clear Time;                                 % Not used in the 5EWA0 labs

%% Calculate

% TODO: insert equations to calculate powers and efficiency

% calculate omega_m continuous
omega_m = n_filtered * 2*pi;

% continuous tfw_DCM
Tf = 0.8 * ((0.001638 * omega_m) + 0.4942);

Td = k_phi*Ia_DCM_filtered;
Ts = Td-Tf;
Ps = Ts .* omega_m;

% combined eff of motor/gen use for loop - for overview
eff_motor = (Ts > 0).*(Td > 0).*Ps./((Va_DCM_filtered.*Ia_DCM_filtered)+(If.*Vf));
% eff_motor(eff_motor == 0) = NaN;

eff_gen = (Ts < 0).*(Td < 0).*(abs(Va_DCM_filtered.*Ia_DCM_filtered))./(abs(Ps)+(If.*Vf));
% eff_gen(eff_gen == 0) = NaN;

eff = eff_motor + eff_gen;

%% Plot data

% TODO: plot the torque vs speed. Think about the axis.
figure;
plot(Ts, omega_m);
title('Torque vs Speed');
xlabel('Speed (rad/s)');
ylabel('Shaft Torque (Nm)');

% TODO: plot the efficiency plots. Use only the datapoints which are
% physically posible. If you calculate an efficiency less than 0% or
% larger than 100% something is wrong.
figure;
hold on;
plot(Ts, eff, 'o');

% plot(n_filtered, eff_gen, 'o');
title('Efficiency vs Speed');
xlabel('Torque (Nm)');
ylabel('Efficiency');
% axis([min(n_filtered) max(n_filtered) 0 1]); % Ensuring efficiency is between 0 and 1


hold off;
% Plot tips:
% - Only plot same type of data in one plot (for example: only powers, no
% combination of power and torque)
% - Try to avoid hold on statements, try to plot in one statement:
% (plot(x,y1,x,y2)
% - Use subplots if you want to see different datatypes on one screen









