%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         
%   Script name:        Lab1_5EWA0_c2_s4_template.m                                                                                    
%   Comments:           Template script for Lab 1 question 2.4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

% TODO: save this script in your own work folder. Remove template from the
% name and add your own name.

%% Measurement data

% TODO: insert measurement data and parameters
n = 40;                               % rotational speed [rev/s]
omega_m = 2*pi*n;                   % rotational speed [rad/s]

If = 0:0.08:0.48;                            % field current setpoints
Va = [2.9, 48.5, 103.3, 154.4, 193.8, 223.3, 244.5];                            % measured armature voltage


%% Calculations 

% TODO: add equation to calculate K_phi
K_phi = Va/omega_m;


%% Plot results

% TODO: add measured data to plot
figure();                           % or use next line for large figure
% figure('units','normalized','outerposition',[0 0 1 1]);                           
plot(If, K_phi,'o');                         % add variables to plot. 'o' symbol is added to represent the data on plot using circles instead of line
grid on;
xlabel('Field current I_f [A]');
ylabel('K_phi');        % add y label

% TODO make a new plot with calculated data 


