clear
clc
close all
addpath('function/')
addpath('plotting function/')
addpath('Ferrofluid interpolation/')
%% Main Code
tic

%% Parameters 
% Interpolation data for EMG water-based series ferrofluid (Ferrotec Inc.)
% URL: https://ferrofluid.ferrotec.com/products/ferrofluid-emg/water/
% EMG   304, 308 408 507  508  509 605  607  700   705  707  708

Chi0_data = [5.03 0.5  0.5  1.63 0.88 0.5  3.02 1.63 12.57 4.04 1.51 0.63];     % Initial Susceptibility
Ms_data   = [27.5 6.6  6.6  11   6.6  3.3  22   11   35.5  22   11   6.6]* 795.7747;      % Saturation Magnetization [A/m]
phi_data  = [4.5  1.2  1.2  2    1.2  0.6  3.9  2    5.8   3.9  2    1.2];      % Magnetic Particle concentration [%]
rho_data  = [1240 1060 1070 1120 1070 1030 1180 1100 1290  1190 1100 1080];     % Density @ 25 C [kg/m^3]

[phi_Chi0_cfit, gof1] = createFit_phi_chi0(Chi0_data, phi_data);
[rho_Chi0_cfit, gof2] = createFit_chi0_rho(Chi0_data, rho_data);
[Ms_phi_cfit, gof3] = createFit_Ms_phi(phi_data, Ms_data);

% Use fitting curves to find Ms from chi0 [User defined value]
chi0    = 0.1;                                                                  % Initial Magnetic Susceptibility of Ferrofluid (User Defined)
rho     = rho_Chi0_cfit(chi0);                                                  % Density Found corresponding to chi0 [kg/m^3]
phi     = phi_Chi0_cfit(chi0); 
Ms      = Ms_phi_cfit(phi);                                                     % Saturation Magnetization Found corresponding to chi0 [A/m]
% Ms      = Ms_rho_cfit(rho)/(1000*4*pi*1e-7)
g       = 9.81 %1.62%9.81%9.806%1e-6%9.806;                                     % Gravitational Acceleration [m/s^2]
mu0     = 4*pi*1e-7;                                                            % Vaccum Pemeability [N/A^2]
gamma   = 3*chi0/Ms; 
sigma   = 0.025;                                                                % Interfacial tension [N/m] water based

%% Execution
parameters  = [g;mu0;rho;Ms;gamma;sigma];

%% User input
mirror_max_length    = 0.275; 
N_count              = 91;
step                 = 0.0001;
M_magnetization      = 1.48;%1.48;                                                   %N52 magnetization, in Tesla
target_smooth_length = 0.1;                                                          %the desired minimal length of smooth surface

x =[0.079013466527307,0.001048291429140,0.023752220839823,0.002937957192350,0.001009665482669,0.009703139524191,0.027966648680831,0.025321322887420];
% X                       = x(1);              %perturbed length
% Z                       = x(2);              %uniform z distance
% lambda                  = x(3);              %length of each repeating pattern
% d                       = x(4);              %uniform vertical spacing
% z1                      = x(5);              %magnet 1 vertical spacing
% z2                      = x(6);              %magnet 2 vertical spacing
% z3                      = x(7);              %magnet 3 vertical spacing
% z4                      = x(8);              %magnet 4 vertical spacing

%% Plotting and outputting values
theta = linspace(-10, 10, 201);

[w,h,mag_count, mag_tot_length, results,mags] = result_plots(x, parameters,N_count,step,M_magnetization,target_smooth_length, theta) ;                                                                         %This function was created to localize the plots somewhere other than the main code.
w
h
mag_count
mag_tot_length

toc