clear
clc
close all
addpath('function')
addpath('plotting function')
addpath('Ferrofluid interpolation')
%% Note
%Rosenweig Instability constraints have not been implemented (@TY)
%implement potential plotting functions (@TY)
%play with penaltyfactor (@eric)
% profile on

%% Main Code
tic
%intcon = 1; We don't need any of our variable be integer
rng default % For reproducibility

A   = [];
b   = [];
Aeq = [];
beq = [];

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

%% Old parameters
% g           = 9.81;                      % Gravitational Acceleration [m/s^2]                                                [parameter(1)]
% % g           = 9.81;
% mu0         = 4*pi*1e-7;                  % Vaccum Pemeability [N/A^2]                                                        [parameter(2)]
% rho         = 1290;                       % Ferrofluid material density [kg/m^3]                                              [parameter(3)]
% % Md    = 28250;                          % Saturation moment of the bulk Ferrofluid [A/m]
% % d     = 10e-9;                      % Nominal Particle diameter of ferrofluid [m]
% % phi   = 1;                          % Volume fraction of Ferrofluid
% % k     = 1.380649e-23;               % Boltzman Constant [m^2*kg/s^2/K]
% % T     = 293;                        % Temperature [K]
% Ms          = 70000 %28250;                      % Saturation magnetization of Ferrofluid mateirial [A/m]                            [parameter(4)]
% chi0        = 0.1;  %0.1 org                      % Initial Magenetic susceptibility of Ferrofluid (EMG 700 with dilution, SI unit)
% gamma       = 3*chi0/Ms;                  %                                                                                   [parameter(5)]
% sigma       = 0.01;                       % Interfacial tension [N/m]

%% Execution
parameters  = [g;mu0;rho;Ms;gamma;sigma];
parameters  = [g;mu0;rho;Ms;gamma;sigma];

%% User input
mirror_max_length    = 0.275; 
N_count              = 91;
step                 = 0.0001;
Cycle_length         = 1;                                                            %pattern frequency (4 blocks of magnet per cycle length)
M_magnetization      = 1.48;%1.48;                                                   %N52 magnetization, in Tesla
target_smooth_length = 0.1;                                                          %the desired minimal length of smooth surface
desired_g            = 9;%10;

Mag_height_lower_lim    = 0.003; %0.003175                                           %m, can be more relaxed depending on the available magnets supplies
Mag_height_upper_lim    = 0.03; %0.0254                                              %m, can be more relaxed depending on the available magnets supplies
Mag_width_lower_lim     = 0.003;%0.003175                                            %m, can be more relaxed depending on the available magnets supplies
Mag_width_upper_lim     = 0.03; %0.0254                                              %m, can be more relaxed depending on the available magnets supplies
Vert_spacing_lower_lim  = 0.001;                                                     %m, can be more relaxed depending on the manufacturing accuracy
Vert_spacing_upper_lim  = Mag_height_upper_lim;                                      %m, can be higher if have taller magnets
Hori_spacing_lower_lim  = 0.001;                                                     %m, can be more relaxed depending on the manufacturing accuracy
Hori_spacing_upper_lim  = Mag_width_upper_lim;                                       %m, can be higher if have girth-ier magnets

%% Setting up lb and ub
%Declare a zero vector of lb and ub
lb                      = zeros(1, 8);
ub                      = zeros(1, 8);

%organization of x will be as follows: x = [uniform height, uniform width, uniform horizontal spacing, perturbed length, vertical spacing between each magnet and mirror surface]
% lb(1)                   = Mag_height_lower_lim;                            %uniform height LOWER bound                                  
% ub(1)                   = Mag_height_upper_lim;                            %uniform height UPPER bound 
% lb(2)                   = Mag_width_lower_lim;                             %uniform width LOWER bound
% ub(2)                   = Mag_width_upper_lim;                             %uniform width UPPER bound
% lb(3)                   = Hori_spacing_lower_lim;                          %uniform horizontal spacing LOWER bound
% ub(3)                   = Hori_spacing_upper_lim;                          %uniform horizontal spacing UPPER bound
% lb(4)                   = 0;                                               %perturbed length LOWER bound
% ub(4)                   = (mirror_length - target_smooth_length)/2 ;       %perturbed length UPPER bound

%the new x = [X, Z, lambda, d, z1, z2, z3, z4]
lb(1)                   = 0;                                                         %perturbed length LOWER bound                                  
ub(1)                   = (mirror_max_length - target_smooth_length)/2;              %perturbed length UPPER bound 

lb(2)                   = Vert_spacing_lower_lim;                                    %Z                             
ub(2)                   = Vert_spacing_upper_lim;                                    %Z                              

lb(3)                   = 4*Hori_spacing_lower_lim + 4*Mag_width_lower_lim ;         %lambda length LOWER bound                       
ub(3)                   = 4*Hori_spacing_upper_lim + 4*Mag_width_upper_lim ;         %lambda length UPPER bound                           

lb(4)                   = Hori_spacing_lower_lim;                                    %uniform horizontal spacing LOWER bound
ub(4)                   = Hori_spacing_upper_lim;                                    %uniform horizontal spacing UPPER bound

lb(5:8)                 = Vert_spacing_lower_lim;                                    %z1-z4
ub(5:8)                 = Vert_spacing_upper_lim;                                    %z1-z4


%% Options settings for G.A.
opts = optimoptions('ga','UseParallel',true);                   %enable "UseParallel" to enable parallet computing and result in ~30% increase
opts.PopulationSize         = 800;
% opts.PopulationSize         = 800*2;                            %200 by default, but 801 is a well proven empirical number
% opts.PopulationSize         = 25601*4;
popSize                     = opts.PopulationSize;
opts.FunctionTolerance      = 1e-15;
opts.ConstraintTolerance    = 1e-12;                            %1e-6 by default
% opts.FitnessLimit          = %[-inf] by default 
opts.PenaltyFactor          = 2; %100 by default              %worth playing with
% opts.MigrationInterval     = 50 %20 by default
% opts.MaxstallGenerations   = 100 %50 by default
opts.MaxGenerations         = 200*length(lb);
opts.EliteCount             = ceil(0.05*popSize);               %{ceil(0.05*PopulationSize)}by default, worth playing with

fun            = @(y) objective(y, parameters, N_count,step,M_magnetization,target_smooth_length,desired_g);
non_linear_con = @(y) nonlcon(y,target_smooth_length,Mag_width_lower_lim,Mag_width_upper_lim,parameters,N_count,step,M_magnetization,desired_g);
% fun = @objective;
% non_linear_con = @nonlcon;

[x,fval,exitflag,output] = ga(fun,length(lb),A,b,Aeq,beq,lb,ub,non_linear_con,opts) %genetic algorithm optimization



%% Plotting and outputting values
plot_x_range            = 0.5;                                                      %simulation of x, from -x to x
plot_y_range            = 0.5;                                                      %simulation of y, from -y to y
plot_mag_resolution     = 1200;
plot_quiver_resolution  = 20;
% % % x = [6.16169716825821e-07	0.0299999993445916	0.0168406873904386	0.00103617017324607	0.0299999875941347	0.00603624152390104	0.0275731371187743	0.0243825982536424]
[w,h,mag_count, mag_tot_length, max_no_trun_error,avg_no_trun_error, max_trun_error,avg_trun_error,mags] = result_plots(x, parameters,N_count,step,M_magnetization,plot_x_range, plot_y_range, plot_mag_resolution, plot_quiver_resolution,target_smooth_length) ;                                                                         %This function was created to localize the plots somewhere other than the main code.
w
h
mag_count
mag_tot_length
max_no_trun_error
avg_no_trun_error
max_trun_error
avg_trun_error
X                       = x(1);              %perturbed length
Z                       = x(2);              %uniform z distance
lambda                  = x(3);              %length of each repeating pattern
d                       = x(4);              %uniform vertical spacing
z1                      = x(5);              %magnet 1 vertical spacing
z2                      = x(6);              %magnet 2 vertical spacing
z3                      = x(7);              %magnet 3 vertical spacing
z4                      = x(8);              %magnet 4 vertical spacing

%% Verifying the magnetic force > 10g constraint
[f_m_mag,violating_index,B_mag] = force_constraint_verify(mag_tot_length,X,N_count,step,parameters,mags,desired_g)
%utilizing fm = mu0*M*delH
%compare and see if H0 same as Hnorm
% [xx, yy, n_x, n_y]          = surface_params(mag_tot_length - 2*X, N_count, step);
% [~,M,Hx,Hy,~,H_norm]        = critical_M(xx,yy,parameters,n_x,n_y,mags);
% h_dir                       = [Hx;Hy]./H_norm;
% Mx                          = h_dir(1,:).*M;
% My                          = h_dir(2,:).*M;
% 
% %Calculating dHdx and dHdy
% % [B_r_step_up, B_z_step_up,~]     = multiB(xx + step*n_x,yy+step*n_y, mags)
% % [B_r_step_down, B_z_step_down,~] = multiB(xx - step*n_x,yy-step*n_y, mags)
% [~,~,B_mag_step_up_wrt_x]   = multiB(xx + step,yy, mags);
% [~,~,B_mag_step_down_wrt_x] = multiB(xx - step,yy, mags);
% [~,~,B_mag_step_up_wrt_y]   = multiB(xx, yy+step, mags);
% [~,~,B_mag_step_down_wrt_y] = multiB(xx, yy-step, mags);
% 
% H_step_up_wrt_r      = B_mag_step_up_wrt_x/mu0;
% H_step_down_wrt_r    = B_mag_step_down_wrt_x/mu0;
% H_step_up_wrt_z      = B_mag_step_up_wrt_y/mu0;
% H_step_down_wrt_z    = B_mag_step_down_wrt_y/mu0;
% % 
% dHdx = 1/(2*step) * (H_step_up_wrt_r - H_step_down_wrt_r);       %equal to zero
% dHdy = 1/(2*step) * (H_step_up_wrt_z - H_step_down_wrt_z);
% 
% f_m = mu0.* (Mx .* dHdx + My .*dHdy)
% if min(abs(f_m))<10*9.81
%     disp('10g constraint violated')
% else
%     disp('constraints satisfied')
% end
% f_m_r          = mu0 *h_dir1(1,:).*M.*[dHdr + dHdz]           %M is calulated from Langevin curve in critical M
% f_m_z          = mu0 *h_dir1(2,:).*M.*[dHdr + dHdz]

toc
% profile viewer
% profsave



