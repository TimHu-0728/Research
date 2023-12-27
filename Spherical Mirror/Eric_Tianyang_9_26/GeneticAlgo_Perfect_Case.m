clear
clc
close all


% openExample('globaloptim/OptimizeANonsmoothFunctionUsingGaExample')

% xi = linspace(-6,2,300);
% yi = linspace(-4,4,300);
% [X,Y] = meshgrid(xi,yi);
% Z = ps_example([X(:),Y(:)]);
% Z = reshape(Z,size(X));
% hold on
% surf(X,Y,Z,'MeshStyle','none')
% colormap 'jet'
% view(-26,43)
% xlabel('x(1)')
% ylabel('x(2)')
% title('ps\_example(x)')

% hold off

%parpool
tic
%intcon = 1; We don't need any of our variable be integer
rng default % For reproducibility
fun = @objective;
A = [];
b = [];
Aeq = [];
beq = [];

r           = 30;                                                               %radius of the mirror's curvature
r_flat      = 15;                                                               %radius of the mirror's projection onto a flat surface
g           = 9.81;                                                             %gravitational constant
N_count     = 91;

xx          = linspace(0, r_flat, N_count);
yy          = surf_func(xx);
r_upper_lim = 25;                                                               %upper limit for radius in the optimization, in meters
r_lower_lim = 0;                                                                %lower limit for radius  

z_upper_lim = 4.5;                                                                %upper limit for height (measured from mirror center), in meters
z_lower_lim = -6; %-(1e-10);                                                               %lower limit for height

I_upper_lim = 18;                                                               %upper limit for ln(current) (measured in ln(Amperes))
I_lower_lim = 11;                                                               %lower limit for ln(current)

w_upper_lim = 2*sqrt(g/r);                                                      %upper limit for rotation rate (measured in radians/s)
w_lower_lim = 0*sqrt(g/r);                                                      %lower limit for rotation rate

IN          = 4;                                                                %number of coils

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

% nonlcon     = [];                                                               %@ellipsecons;
opts = optimoptions('ga','UseParallel',false);
opts.PopulationSize         = 802; %200 by default, 802 is optimal (perfect case)
popSize                     = opts.PopulationSize;
opts.FunctionTolerance      = 1e-15;
opts.ConstraintTolerance    = 1e-12;
 %1e-6 by default
% opts.FitnessLimit          = %[-inf] by default 
opts.PenaltyFactor         = 2; %100 by default
% opts.MigrationInterval      = 50 %20 by default
% opts.MaxstallGenerations    = 100 %50 by default
opts.MaxGenerations         = 200*length(lb);
opts.EliteCount             = ceil(0.05*popSize); %{ceil(0.05*PopulationSize)}by default

[x,fval,exitflag,output] = ga(fun,IN * 3 + 1,A,b,Aeq,beq,lb,ub,@nonlcon, opts)         %genetic algorithm optimization

%initializing variables for plotting
% IN              = (max(size(x)) - 1) / 3;                                                             %number of coils
rz              = zeros(IN, 2);                                                                       %initialzing vectors for massForcePotential inputs - [r, z]
I               = zeros(IN, 1);                                                                       %I in massForcePotential
sum_RI          = 0;                                                                                  %initializing the sum_RI measure to check
for i = 1:IN                                                                                          %generating massForcePotential inputs
    I(i, 1)     = exp(x(i * 3 - 1));
    rz(i, 2)    = x(i * 3);
    rz(i, 1)    = x(i * 3 + 1);
    x(3 * i - 1)= exp(x(i * 3 - 1));
    sum_RI      = sum_RI + abs(x(3 * i - 1)*x(3 * i + 1));
end

%Things to check
% disp(x)                                                                                               %Displaying the chosen optimizied values
sum_RI          = sum_RI                                                                              %Limit is +-1500000, I*r proportional to mass of the superconductor
% weighted_PI_difference    = fval                                                                      %user defined objective function value after subbing in x

%Plotting


result_plots(x, r_flat, I, rz, N_count);                                                                             %This function was created to localize the plots somewhere other than the main code.
toc


function [c ceq] = nonlcon(y)
    % Parameters using EMG 700 from Ferrotec corp. Website
    % URL: https://ferrofluid.ferrotec.com/products/ferrofluid-emg/water/emg-700-sp/
    g       = 9.806;                      % Gravitational Acceleration [m/s^2]
    mu0     = 4*pi*1e-7;                  % Vaccum Pemeability [N/A^2]
    rho     = 1290;                       % Ferrofluid material density [kg/m^3]
    Ms      = 28250;                      % Saturation magnetization of Ferrofluid mateirial [A/m]          
    chi0    = 0.1;                        % Initial Magenetic susceptibility of Ferrofluid (EMG 700 with dilution, SI unit)
    gamma   = 3*chi0/Ms; 
    sigma   = 0.01;                       % Interfacial tension [N/m]
    r_flat  = 15;
    N_count = 901;
    step    = 0.0001;
    [r, z, nr, nz] = surface_params(r_flat, N_count, step);

    parameters = [g;mu0;rho;Ms;gamma;sigma];
    I          = [y(2);y(5);y(8);y(11)];
    co         = [y(3) y(4);y(6) y(7);y(9) y(10);y(12) y(13)];
    
    % c(1) = y(3) - surf_fun(y(4))
    % IN  = 4;
    %     for i = 1:IN
    %     c(i) = y(3*i) - surf_func(y(3*i +1));
    %     end
    % ceq = [];
    
    [Mc,M]   = Critical_M(r',z',nr,nz,I,co,parameters);
    

    IN  = 4;
    c   = zeros(IN+N_count);
        for i = 1:IN
            c(i) = y(3*i) - (30 - sqrt(900 - y(3*i +1) .^ 2));
        end

        for i = 1:N_count
            c(4+i) = M(i)-Mc(i);
        end

    ceq = [];
end

function obj = objective(y)
    IN = (max(size(y)) - 1) / 3;                                                                       %number of coils
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
    top          = massForcePotential(y(1),xx,yy,I,rz) ...
        - massForcePotential(y(1),X_ideal,Y_ideal,I,rz);
    step         = 0.0001;
    bottom       = 1/(2*step) * (massForcePotential(y(1),xx+step*nx,yy+step*ny,I,rz) ...
        - massForcePotential(y(1),xx-step*nx,yy-step*ny,I,rz));
    error        = top./bottom;
    L2norm       = sqrt(sum(sum(error.^2)));                                                           %L2norm of the error (magnitude)
    Linf         = max(max(abs(error)));                                                               %to account for sudden spikes
    gamma        = 1;                                                                                  %tuning parameter for Linf term
    obj          = L2norm/sqrt(max(size(error))) + gamma *Linf;
%     obj          = L2norm + gamma *Linf;
end
