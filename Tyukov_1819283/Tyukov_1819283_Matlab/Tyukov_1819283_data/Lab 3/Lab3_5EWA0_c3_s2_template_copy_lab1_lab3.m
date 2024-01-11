clear all;
close all;
clc;

% TODO: save this script in your own work folder. Remove template from the
% name and add your own name.
% COPY OF LAB 1

%% Measurement data

% TODO: fill measurement results in data vectors
n = [1,5,10,20,30,40];                             % rotational speed [rev/s]
omega_m = 2*pi*n;                         % rotational speed [rad/s]
Ia = [0.44, 0.59, 0.67, 0.80, 0.91, 0.95];                            % measured armature current
k_phi = 0.968;                           % calculated value from previous section


%% Calculations 

% TODO: add equation to calculate friction torque T_f and mechanical power
% delivered by the DC machine.
T_f = k_phi*Ia;

%% Plot results

figure('units','normalized','outerposition',[0 0 1 1]);

% TODO: add data to plot the speed vs the friction torque
plot(omega_m, T_f,'o');                         % add variables to plot. 'o' symbol is added to represent the data on plot using circles instead of line
grid on;
xlabel('add clear x-label');        % add x label
ylabel('add clear y-label');        % add y label

%% Data fit 
% TODO: Add fit to measured curve of friction torque Tfw using the fit function.
%
% fitobject = fit(x,y,fitType) creates the fit to the data in x and y with the model specified by fitType.
% fitType is a model type to fit, which can be chosen from, for example:
%
% 1. 'poly1' - Linear polynomial curve
% 2. 'poly2' - Quadratic polynomial curve
% 3. 'poly3' - Cubic polynomial curve
%
% Think about the order of polynomial order based on the type of dependency you have. 

% 2. plot fitting result. 
%
% Note: Instead of line vectors, the x and y axis vectors should be transformed 
% in columnn vectors : x' and y'.

f = fit(omega_m', T_f','poly1')

% Plot fit result in same window as the measurement data. 
hold all;
plot(f);
xlabel('omega_m');        % add x label
ylabel('T_f');        % add y label
%% Plot results

% TODO: Finalise plots with legend, axis labels, units etc.
