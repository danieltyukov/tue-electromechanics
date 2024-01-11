format long
%% Parameters
Va = 600;
Ia = 50;
Ra = 0.5;
Ra2 = 5;

Vf = 100;
Rf = 100;
Rf2 = 100;

revs = 1500;
Tfw = 20;


%% Question 1
E = Va-Ra*Ia;
w = (revs/60)*2*pi;
kp=E/w;
Td = kp*Ia;

%% Question 2
Tload = Td-Tfw;

%% Question 3
Vf2 = 120;
If1 = Vf/Rf;
If2 = Vf2/Rf;
kp2 = (kp*If2)/If1;

Ia2 = Td/kp2;
E2 = kp2*w;
Va2 = E2+Ra*Ia2;

%% Question 4
Va_shunt = 20;
Va_shunt2 = 40;
Rf_shunt = 100;
Ra_shunt = 5;

If_shunt = Va_shunt/Rf_shunt;
Ia_shunt = Va_shunt/Ra_shunt;
Ia_shunt2 = Va_shunt2/Ra_shunt;
Td_shunt1 = 10;
kp_shunt = Td_shunt1/Ia_shunt;

kp_shunt2 = (kp_shunt*Ia_shunt2)/Ia_shunt;
Td_shunt2 = Ia_shunt2*kp_shunt;

%% Question 5
w2 = 400/51;

%% Question 6
Ia_ps = (8-w2) + Va_shunt2/Rf_shunt;

%% Question 7
Va_series = 300;
Ia_series = 60;
Rf_series = 0.8;
Ra_series = 0.3;
Tload_series = 130;
w_series = 2*pi*(1000/60);

Td_series = (Ia_series*(Va_series-(Rf_series+Ra_series)*Ia_series))/(w_series);

Tfw_series = Td_series - Tload_series;


%% Question 8
Pout = Tload_series*w_series;
Pin = Va_series*Ia_series;
eta = 100*Pout/Pin;

%% Question 9
Rs_series = 0.2;
w_series2 = (Ia_series*(Va_series-(Rf_series+Rs_series+Ra_series)*Ia_series))/(Td_series);


%% Question 10
Rl_gen = 3;
Tin_gen = 10;
kp_gen = (Va_series-Ia_series*(Ra_series+Rf_series))/w_series;
w_gen = Ia_series*(Rf_series+Ra_series+Rl_gen)/kp_gen