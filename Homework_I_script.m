format long
%% PROVIDED CONSTANTS
g = 1*10^-3; % mm
hc = 18 * 10^-3; %mm
hl = 30*10^-3; % mm
hm = 15*10^-3; % mm
ht = 10*10^-3; % mm
lm = 40*10^-3; % mm
w1 = 15*10^-3; % mm
w2 = 22.5*10^-3; % mm

wc = 8*10^-3; %% mm
ws = 50*10^-3; % mm
wp = 5*10^-3;  % mm
wt = 87.5*10^-3; % mm

lt1 = 5*10^-3; % mm
lt2 = 68.75*10^-3; % mm
lt3 = 5*10^-3; % mm
D = 10*10^-3;  % mm
SF = 0.95; 
Br = 1.3;      % T
Hc = 1000*10^3; % A/m
N1 = 900;      % Turns
N2 = 625;      % Turns
mu_r = 1000; 
mu_0 = 4*pi*10^-7;

% -< comment mark means verified correct

len1 = wp+(w1/2)+(hm/2)+hl;
len2 = (wp+lm-ws)+(hm/2)+(w2/2);

%% exercise 3a
mu_pr = Br/(Hc*mu_0); % 

%% exercise 4
% magnet
R_pm=lm/(mu_0*mu_pr*hm*D); %

% left side
R_st1 = len1/(mu_0*mu_r*w1*D); %
R_g1  = g/(mu_0*w1*D); %
R_sb1 = lt1/(mu_0*mu_r*w1*D); %

% middle
R_sb2 = lt2/(mu_0*mu_r*ht*D); %

% right side
R_st2 = hl/(mu_0*mu_r*w2*D*SF); % SF
R_st3 = len2/(mu_0*mu_r*hm*D); % 
R_g2  = g/(mu_0*w2*D); %
R_sb3 = lt3/(mu_0*mu_r*w2*D); %

% x=R_sb1+R_sb2+R_sb3;

%% exercise 5
R_eq = R_pm + R_st1 + R_st2 + R_st3 + R_g1 + R_g2 + R_sb1 + R_sb2 + R_sb3;

flux = (Hc*lm)/R_eq;

B_P1 = flux/(w1*D); % ECF
B_P2 = flux/(hm*D); %
B_P3 = flux/(w2*D*SF); %

%% exercise 6
L1 = (N1^2)/R_eq; % sus
L2 = (N2^2)/R_eq;
M = (N1*N2)/R_eq;

%% exercise 7
mmf = Hc*lm;

syms D_g;
Peq = 1/(D_g/(mu_0*w1*D) + D_g/(mu_0*w2*D) + lm/(mu_0*hm*D));
dP = diff(Peq,D_g);
F_d = 1/2*(mmf)^2*dP;

F=double(subs(F_d, D_g, g));

% Force is upwards, closing the gap

%% exercise 8

mmf_plusten = Hc*lm +N2*10;
mmf_minusten = Hc*lm - N2*10;

F_plus10 = double(subs(0.5*dP*(mmf_plusten)^2, D_g, g));
F_minus10 = double(subs(0.5*dP*(mmf_minusten)^2, D_g, g));

%% exercise 9

Btop = flux/(hm*D);
Htop = Btop*(400+exp(3.9*Btop));
mmf_top = Htop*(w1/2 + wp+ w2/2);

Bleft = flux/(w1*D);
Hleft = Bleft*(400+exp(3.9*Bleft));
mmf_left = Hleft*(hl+hm/2);

Bbottom= flux/(ht*D);
Hbottom = Bbottom*(400+exp(3.9*Bbottom));
mmf_bottom = Hbottom*wt;

Bright = flux/(w2*D*SF);
Hright = Bright*(400+exp(3.9*Bright));
mmf_right = Hright*(hl+hm/2);

Bair1 = 1.8;
mmf_air1 = g*(Bair1/mu_0);

Bair2 = flux/(w2*D);
mmf_air2 = g*(Bair2/mu_0);

mmf_pm = R_pm*flux;

mmf_sum = mmf_top+mmf_left+mmf_bottom+mmf_right+mmf_air1+mmf_air2+mmf_pm;

I1 = ((Hc*lm)-mmf_sum)/N1;
%% exercise 10

