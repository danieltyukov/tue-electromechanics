format long; %increases precision

%% Parameters DONE
V_Lbr=20;
I_Lbr=572.2;
P_tbr=10000;
V_Ldc=1;
I_Ldc=83.3;
f_br=40;
f=160;


%% Exercise 2 DONE
R_eqbr = P_tbr/(3*I_Lbr^2);
Z_br = V_Lbr/(I_Lbr*sqrt(3));
X_br = sqrt((Z_br)^2-(R_eqbr)^2);
R1 = 1.05*(V_Ldc/(2*I_Ldc));
R2 = R_eqbr-R1;
X_eq = (X_br*f)/f_br;
x1 = X_eq/2;
X1_2_norm = X_eq/2;

%% Exercise 4 DONE
s=2.4279*10^(-4);
V_nl=400/sqrt(3);
I_lnl=234.638;
P_tnl=(19.346*10^3)/3;
pf=P_tnl/(V_nl*I_lnl);
theta_nl=-acos(pf);

I_lnl_ph=I_lnl*exp(1i*theta_nl);
E1=V_nl-I_lnl_ph*(R1+1i*x1);
I_2_lnl_ph=E1/((R2/s)+1i*x1);

Rc=((abs(E1))^2)/(V_nl*I_lnl*cos(theta_nl)-(((I_lnl)^2)*R1)-((R2*((abs(I_2_lnl_ph)))^2)/s));
Xm=((abs(E1))^2)/(V_nl*I_lnl*sin(-theta_nl)-(((I_lnl)^2)*X1_2_norm)-((abs(I_2_lnl_ph))^2)*X1_2_norm);

%% Exercise 5 DONE
Pd_nl=3*((abs(I_2_lnl_ph))^2)*((1-s)/s)*R2;
Pfw_nl=Pd_nl;

%% Exercise 6 DONE
s_1 = (4800-4757.6)/4800;
Z_prl = (Rc*1i*Xm)/(Rc+1i*Xm);
Z_ser=R2/s_1+1i*x1;
Zbruh = 1/(1/Z_prl +1/Z_ser);
rad2deg(angle(Z_prl));
rad2deg(angle(Z_ser));
Z_fin = Zbruh + R1 +1i*x1;
I_input=V_nl/Z_fin;
I_input_theta_rad=angle(I_input);
I_input_theta_deg=rad2deg(angle(I_input));
I_input_mag=abs(I_input);

%% Exercise 7 DONE
PF_norm=cos(I_input_theta_rad);

%% Exercise 8 DONE
P_cu_s=3*((I_input_mag)^2)*R1;
P_in=3*V_nl*I_input_mag*PF_norm;
Pd=310*1000+Pfw_nl;
Pg=Pd/(1-s_1);
P_core=P_in-Pg-P_cu_s;
P_cu_r=Pg*s_1;

%% Exercise 9
Efficiency=((Pd-Pfw_nl)*100)/P_in;