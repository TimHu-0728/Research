clc
clear
close all
format long

%% 6 coil
optimal_soln = [0.591722704	2467641.733	3.188185244	16.05227601	8908374.854	-5.999995108	6.779684702	20735.27118	-5.99999902	24.81730996	4051.644025	-5.999994109	24.81730709	11.7910758	-5.999999115	25	2.916066099	-5.999999969	25];
parameters = [9.806 1.25663706143592e-06 1045.639438 2849.727161 0.1 0.025];

x0 = [optimal_soln parameters];
f0 = errorFcn_6coils(x0);

% sa
dx =0.1*ones(size(x0));
x_sen=zeros(length(x0),1);
parfor i = 1:length(x0) 
    x_sen(i) = (abs((errorFcn_6coils(changem(x0,x0(i)+dx(i),x0(i))))-f0))/dx(i);
end

w_sen=x_sen(1);
I1_sen=x_sen(2);
z1_sen=x_sen(3);
r1_sen=x_sen(4);
I2_sen=x_sen(5);
z2_sen=x_sen(6);
r2_sen=x_sen(7);
I3_sen=x_sen(8);
z3_sen=x_sen(9);
r3_sen=x_sen(10);
I4_sen=x_sen(11);
z4_sen=x_sen(12);
r4_sen=x_sen(13);
I5_sen=x_sen(14);
z5_sen=x_sen(15);
r5_sen=x_sen(16);
I6_sen=x_sen(17);
z6_sen=x_sen(18);
r6_sen=x_sen(19);
g_sen=x_sen(20);
rho_sen=x_sen(22);
Ms_sen=x_sen(23);
chi0_sen=x_sen(24);
sigma_sen=x_sen(25);


%% 5 coils

optimal_soln = [0.585459507	6016241.326	-5.908527641	7.004079989	7399466.245	3.361582194	17.29333946	689920.7492	4.429855921	22.26700704	52.42660048	-5.999999958	24.99999999	4.678936735	-5.999999987	25];
parameters = [9.806 1.25663706143592e-06 1045.639438 2849.727161 0.1 0.025];

x0 = [optimal_soln parameters];
f0 = errorFcn_5coils(x0);

% sa
dx =0.1*ones(size(x0));
x_sen=zeros(length(x0),1);
parfor i = 1:length(x0) 
    x_sen(i) = (abs((errorFcn_5coils(changem(x0,x0(i)+dx(i),x0(i))))-f0))/dx(i);
end

w_sen=x_sen(1);
I1_sen=x_sen(2);
z1_sen=x_sen(3);
r1_sen=x_sen(4);
I2_sen=x_sen(5);
z2_sen=x_sen(6);
r2_sen=x_sen(7);
I3_sen=x_sen(8);
z3_sen=x_sen(9);
r3_sen=x_sen(10);
I4_sen=x_sen(11);
z4_sen=x_sen(12);
r4_sen=x_sen(13);
I5_sen=x_sen(14);
z5_sen=x_sen(15);
r5_sen=x_sen(16);
g_sen=x_sen(17);
rho_sen=x_sen(19);
Ms_sen=x_sen(20);
chi0_sen=x_sen(21);
sigma_sen=x_sen(22);

%% 4 coils

optimal_soln = [0.584659474	5640352.185	-5.949607568	6.92109472	7372279.22	3.155655999	17.06239268	6.052319957	-5.999999984	25	6.091067175	-5.999999987	25];
parameters = [9.806 1.25663706143592e-06 1045.639438 2849.727161 0.1 0.025];

x0 = [optimal_soln parameters];
f0 = errorFcn_4coils(x0);

% sa
dx =0.1*ones(size(x0));
x_sen=zeros(length(x0),1);
parfor i = 1:length(x0) 
    x_sen(i) = (abs((errorFcn_4coils(changem(x0,x0(i)+dx(i),x0(i))))-f0))/dx(i);
end

w_sen=x_sen(1);
I1_sen=x_sen(2);
z1_sen=x_sen(3);
r1_sen=x_sen(4);
I2_sen=x_sen(5);
z2_sen=x_sen(6);
r2_sen=x_sen(7);
I3_sen=x_sen(8);
z3_sen=x_sen(9);
r3_sen=x_sen(10);
I4_sen=x_sen(11);
z4_sen=x_sen(12);
r4_sen=x_sen(13);
g_sen=x_sen(14);
rho_sen=x_sen(16);
Ms_sen=x_sen(17);
chi0_sen=x_sen(18);
sigma_sen=x_sen(19);

%% 10 coils

optimal_soln = [0.591536755	2692050.985	-4.469972142	4.998345065	3960124.114	3.118144946	16.57518241	7296173.487	-5.887909122	24.48747852	2858326.745	4.489339481	24.99998944	202576.8289	-5.940944287	24.64098385	21089.39999	4.217669586	23.3367982	6733656.579	-5.864280649	9.58054862	16630.12937	2.672357685	23.37035276	695.5871406	-5.999983534	25	162.3844021	-4.520244486	23.71319643];
parameters = [9.81 1.25663706143592e-06 1045.639438 2849.727161 0.1 0.025];

x0 = [optimal_soln parameters];
f0 = errorFcn_10coils(x0);

% sa
dx =0.0001*abs(x0);
x_sen=zeros(length(x0),1);
parfor i = 1:length(x0) 
    x_sen(i) = (abs((errorFcn_10coils(changem(x0,x0(i)+dx(i),x0(i))))-f0))/dx(i);
end

w_sen=x_sen(1)*1e9;
I1_sen=x_sen(2)*1e6;
z1_sen=x_sen(3)*1e6;
r1_sen=x_sen(4)*1e6;
I2_sen=x_sen(5)*1e6;
z2_sen=x_sen(6)*1e6;
r2_sen=x_sen(7)*1e6;
I3_sen=x_sen(8)*1e6;
z3_sen=x_sen(9)*1e6;
r3_sen=x_sen(10)*1e6;
I4_sen=x_sen(11)*1e6;
z4_sen=x_sen(12)*1e6;
r4_sen=x_sen(13)*1e6;
I5_sen=x_sen(14)*1e6;
z5_sen=x_sen(15)*1e6;
r5_sen=x_sen(16)*1e6;
I6_sen=x_sen(17)*1e6;
z6_sen=x_sen(18)*1e6;
r6_sen=x_sen(19)*1e6;
I7_sen=x_sen(20)*1e6;
z7_sen=x_sen(21)*1e6;
r7_sen=x_sen(22)*1e6;
I8_sen=x_sen(23)*1e6;
z8_sen=x_sen(24)*1e6;
r8_sen=x_sen(25)*1e6;
I9_sen=x_sen(26)*1e6;
z9_sen=x_sen(27)*1e6;
r9_sen=x_sen(28)*1e6;
I10_sen=x_sen(29)*1e6;
z10_sen=x_sen(30)*1e6;
r10_sen=x_sen(31)*1e6;
g_sen=x_sen(32)*1e9;
rho_sen=x_sen(34)*1e9;
Ms_sen=x_sen(35)*1e9;
chi0_sen=x_sen(36)*1e9;
sigma_sen=x_sen(37)*1e9;