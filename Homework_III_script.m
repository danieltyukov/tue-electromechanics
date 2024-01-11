format long
%% Parameters
f=60;
V=400;
Van = V/sqrt(3);
rpm=900;
wm = 900*2*pi/60;
we = 2*pi*f;

%% Exercise 1
p = 120*f/rpm;
pp = p/2;

%% Exercise 3
Pd = 2500;
%cos(theta) = 0.9735;
PF = 0.9851;
theta = -acos(PF);

Ia = (Pd*sqrt(3))/(3*V*PF);

%% Exercise 5
Rab = 8.72;
Rac = 8.82;
Rbc = 8.78;

Lab = 0.07061;
Lac = 0.08146;
Lbc = 0.07811;

R_avg = (Rab+Rac+Rbc)/6;
L_avg = (Lab+Lac+Lbc)/6;

%% Exercise 6
Ef = Van-Ia*exp(1i*theta)*(1i*we*L_avg+R_avg);
mag_Ef = abs(Ef);

%% Exercise 7
delta = angle(Ef);

%% Exercise 8
beta = delta - theta;
Td = (3*mag_Ef*Ia*cos(beta))/wm;

%% Exercise 9
P_cu = 3*R_avg*Ia^2
P_CFW = 0.05*wm*Td
Pout = 0.95*Td*wm;
Pin = 2500;
eta = 100*Pout/Pin