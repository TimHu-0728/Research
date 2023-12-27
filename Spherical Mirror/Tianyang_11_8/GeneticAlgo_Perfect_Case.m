clear
clc
close all


%% 
% Interpolation data for EMG water-based series ferrofluid (Ferrotec Inc.)
% URL: https://ferrofluid.ferrotec.com/products/ferrofluid-emg/water/
% EMG   304, 308 408 507  508  509 605  607  700   705  707  708
Chi0_data = [5.03 0.5  0.5  1.63 0.88 0.5  3.02 1.63 12.57 4.04 1.51 0.63];     % Initial Susceptibility
Ms_data   = [27.5 6.6  6.6  11   6.6  3.3  22   11   35.5  22   11   6.6] * 795.7747;      % Saturation Magnetization [A/m]
phi_data  = [4.5  1.2  1.2  2    1.2  0.6  3.9  2    5.8   3.9  2    1.2];      % Magnetic Particle concentration [%]
rho_data  = [1240 1060 1070 1120 1070 1030 1180 1100 1290  1190 1100 1080];     % Density @ 25 C [kg/m^3]


% Find linear fit curves

[phi_Chi0_cfit, gof1] = createFit_phi_chi0(Chi0_data, phi_data);
[rho_Chi0_cfit, gof2] = createFit_chi0_rho(Chi0_data, rho_data);
[Ms_phi_cfit, gof3] = createFit_Ms_phi(phi_data, Ms_data);

% Use fitting curves to find Ms from chi0 [User defined value]
chi0    = 0.1;                                                                  % Initial Magnetic Susceptibility of Ferrofluid (User Defined)
rho     = rho_Chi0_cfit(chi0);                                                  % Density Found corresponding to chi0 [kg/m^3]
phi     = phi_Chi0_cfit(chi0); 
Ms      = Ms_phi_cfit(phi);                          % Saturation Magnetization Found corresponding to chi0 [A/m]
% Ms      = Ms_rho_cfit(rho)/(1000*4*pi*1e-7)
g       = 9.806;                                                             % Gravitational Acceleration [m/s^2]
mu0     = 4*pi*1e-7;                                                         % Vaccum Pemeability [N/A^2]
gamma   = 3*chi0/Ms; 
sigma   = 0.025;                                                             % Interfacial tension [N/m]

parameters = [g;mu0;rho;Ms;gamma;sigma];

%parpool
tic
%intcon = 1; We don't need any of our variable be integer
rng default % For reproducibility
fun = @(y) objective(y,parameters);
nonl_constraint = @(y) nonlcon(y,parameters);

A = [];
b = [];
Aeq = [];
beq = [];  

r           = 30;                                                               %radius of the mirror's curvature
r_flat      = 15;                                                               %radius of the mirror's projection onto a flat surface
% g           = 9.81;                                                             %gravitational constant
N_count     = 91;

xx          = linspace(0, r_flat, N_count);
yy          = surf_func(xx);
r_upper_lim = 25;                                                               %upper limit for radius in the optimization, in meters
r_lower_lim = 0;                                                                %lower limit for radius  

z_upper_lim = 4.5;                                                                %upper limit for height (measured from mirror center), in meters
z_lower_lim = -6; %-(1e-10);                                                               %lower limit for height

I_upper_lim = 20;%30;                                                               %upper limit for ln(current) (measured in ln(Amperes))
I_lower_lim = 0;                                                               %lower limit for ln(current)

w_upper_lim = 2*sqrt(g/r);                                                      %upper limit for rotation rate (measured in radians/s)
w_lower_lim = 1*sqrt(g/r);                                                      %lower limit for rotation rate

IN          = 5;                                                                %number of coils

lb                = zeros(1, 3 * IN + 1);                                       %pre-generating lower and upper boundary constraints
ub                = zeros(1, 3 * IN + 1);
lb(1)             = w_lower_lim;                                                %adding rotation rate, which does not repeat with current
ub(1)             = w_upper_lim;
for i = 1:IN                                                                    %generating boundary constraints for genetic algorithm
    lb(3 * i - 1) = I_lower_lim;
    lb(3 * i)     = z_lower_lim;
    lb(3 * i + 1) = r_lower_lim;
    ub(3 * i - 1) = I_upper_lim;
    ub(3 * i)     = z_upper_lim;
    ub(3 * i + 1) = r_upper_lim;
end
% 
% x0 = [0.584158342	log(929320.2649)	-5.185174203	6.209628052	log(1206691.227)	2.992480404	16.27979142	log(78420.44437)	-4.188205749	24.99999558	log(59875.95195)	-5.415838695	24.99995462];
% nonlcon     = [];                                                               %@ellipsecons;
opts = optimoptions('ga','UseParallel',true);
% opts = optimoptions('particleswarm','SwarmSize',100,'UseVectorized',true,'HybridFcn',@fmincon,'InitialSwarmMatrix',x0);
opts.PopulationSize         = 20000; %200 by default, 802 is optimal (perfect case)
popSize                     = opts.PopulationSize;
opts.FunctionTolerance      = 1e-15;
opts.ConstraintTolerance    = 1e-13;
% opts.InitialPopulationMatrix = x0;
 %1e-6 by default
% opts.FitnessLimit          = %[-inf] by default 
opts.PenaltyFactor          = 100; %100 by default
% opts.MigrationInterval      = 50 %20 by default
% opts.MaxstallGenerations    = 100 %50 by default
opts.MaxGenerations         = 200*length(lb);
opts.EliteCount             = ceil(0.05*popSize); %{ceil(0.05*PopulationSize)}by default

[x,fval,exitflag,output] = ga(fun,IN * 3 + 1,A,b,Aeq,beq,lb,ub,nonl_constraint, opts)         %genetic algorithm optimization

%initializing variables for plotting
% IN              = (max(size(x)) - 1) / 3;                                                             %number of coils
rz              = zeros(IN, 2);                                                                       %initialzing vectors for massForcePotential inputs - [r, z]
I               = zeros(IN, 1);                                                                       %I in massForcePotential
sum_RI          = 0;                                                                                  %initializing the sum_RI measure to check
w               = x(1);
for i = 1:IN                                                                                          %generating massForcePotential inputs
    I(i, 1)     = exp(x(i * 3 - 1));
    rz(i, 2)    = x(i * 3);
    rz(i, 1)    = x(i * 3 + 1);
    x(3 * i - 1)= exp(x(i * 3 - 1));
    sum_RI      = sum_RI + abs(x(3 * i - 1)*x(3 * i + 1));
end

%Things to check
% disp(x)                                                                                               %Displaying the chosen optimizied values
sum_RI                                                                                   %Limit is +-1500000, I*r proportional to mass of the superconductor
w
I
rz

% weighted_PI_difference    = fval                                                                      %user defined objective function value after subbing in x
toc



%% Plotting
result_plots(x, r_flat, I, rz, N_count,parameters);                                                                             %This function was created to localize the plots somewhere other than the main code.


%% Functions
function [c,ceq] = nonlcon(y,parameters)
    
    r_flat         = 15;
    N_count        = 31;
    step           = 0.0001;
    [r, z, nr, nz] = surface_params(r_flat, N_count, step);
    I              = exp([y(2);y(5);y(8);y(11)]);
    co             = [y(4) y(3);y(7) y(6);y(10) y(9);y(13) y(12)];
    
    % c(1) = y(3) - surf_fun(y(4))
    % IN  = 4;
    %     for i = 1:IN
    %     c(i) = y(3*i) - surf_func(y(3*i +1));
    %     end
    % ceq = [];
    
    % Calculate critical magnetization
    [Mc,M]   = Critical_M(r',z',nr',nz',I,co,parameters);
    

    IN  = 5;
    c   = zeros(IN+N_count,1);
    % Constraints for coil being under the mirror surface
    c(1:IN)=[y(3)-(30 - sqrt(900 - y(4) .^ 2));
             y(6)-(30 - sqrt(900 - y(7) .^ 2))
             y(9)-(30 - sqrt(900 - y(10) .^ 2))
             y(12)-(30 - sqrt(900 - y(13) .^ 2))];

    % Constraints for Material Magnetization smaller than critical magnetization along the surface (Avoid Rosensweig Instability)
    c(IN+1:end) = M - Mc;

    ceq = [];
end

function obj = objective(y,parameters)
    % % Parameters using EMG 700 from Ferrotec corp. Website
    % % URL: https://ferrofluid.ferrotec.com/products/ferrofluid-emg/water/emg-700-sp/
    % g     = 9.806;                      % Gravitational Acceleration [m/s^2]
    % mu0   = 4*pi*1e-7;                  % Vaccum Pemeability [N/A^2]
    % rho   = 1290;                       % Ferrofluid material density [kg/m^3]
    % Ms    = 28250;                      % Saturation magnetization of Ferrofluid mateirial [A/m]          
    % chi0  = 1;                          % Initial Magenetic susceptibility of Ferrofluid (EMG 700 with dilution, SI unit)
    % gamma = 3*chi0/Ms; 
    % 
    % parameters = [g;mu0;rho;Ms;gamma];
    w  = y(1);
    IN = (max(size(y)) - 7) / 3;                                                                       %number of coils
    parameters = parameters(1:end-1);
    rz = zeros(IN, 2);                                                                                 %initialzing vectors for massForcePotential inputs - [r, z]
    I = zeros(IN, 1);                                                                                  %I in massForcePotential
    for i = 1:IN                                                                                       %generating massForcePotential inputs
        I(i, 1) = exp(y(i * 3 - 1));
        rz(i, 2) = y(i * 3);
        rz(i, 1) = y(i * 3 + 1);
    end
    N_count      = 91;
    r_flat       = 15;                                                                                  %radius of the mirror's projection onto a flat surface
    step         = 0.0001;                                                                                     %step size
    [xx, yy, nx, ny] = surface_params(r_flat, N_count, step);                                          %finding points to measure potential, and normal vector away from desired surface at those points
    X_ideal      = zeros(1, N_count);
    Y_ideal      = zeros(1, N_count);
    top          = massForcePotential(w,xx,yy,I,rz,parameters) ...
        - massForcePotential(w,X_ideal,Y_ideal,I,rz,parameters);
    step         = 0.0001;
    bottom       = 1/(2*step) * (massForcePotential(w,xx+step*nx,yy+step*ny,I,rz,parameters) ...
        - massForcePotential(w,xx-step*nx,yy-step*ny,I,rz,parameters));
    error        = top./bottom;
    L2norm       = sqrt(sum(sum(error.^2)));                                                           %L2norm of the error (magnitude)
    Linf         = max(max(abs(error)));                                                               %to account for sudden spikes
    gamma2        = 1;                                                                                  %tuning parameter for Linf term
    obj          = L2norm/sqrt(max(size(error))) + gamma2 *Linf;
%     obj          = L2norm + gamma *Linf;
end



