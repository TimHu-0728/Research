function PI = massForcePotential(w,r,z,I,co,parameters)
% massForcePotential calculates 2D mass force potentials PI, PI_m, PI_w, PI_g, 
% and magenetic flux density Br, Bz in axisymmetry cylindrical coordinates 
% from input: 
%         1. [r,z]: Position Vector investigated             [m]
%         2. w: Angular velocity  of the rotating mirror     [rad/s]
%         3. I: Currents of coil magnets                     [A]
%         4. co: Coordinates of coil magnets                 [m]

%% Parameters
g     = parameters(1);                       % Gravitational Acceleration [m/s^2]
mu0   = parameters(2);                       % Vaccum Pemeability [N/A^2]
rho   = parameters(3);                       % Ferrofluid material density [kg/m^3]                     
Ms    = parameters(4);                       % Saturation Magnetization of ferrofluid [A/m]
gamma = parameters(5); 

%% Calculating Magnetic Field H
[Br,Bz] = multiB(r, z, I, co);       % Magnetic flux density [T]
B       = sqrt(Br.^2+Bz.^2);         % Magnetic flux density norm [T]
H_norm  = 1/mu0*B;                   % Magnetic field [A/m] (M = [0;0] assumed)

%% Calculate Mass Force Potentials PI
% Langevin formula
PI_m  = -mu0/rho*Ms*(lnsinh(H_norm * gamma)/gamma);
PI_g  = g*z;
PI_w  = -w^2*r.^2/2;
PI    = PI_g+PI_w+PI_m;
end



% function [PI,PI_m,PI_w,PI_g,Br,Bz] = massForcePotential(w,r,z,I,co,parameters)
% % massForcePotential calculates 2D mass force potentials PI, PI_m, PI_w, PI_g, 
% % and magenetic flux density Br, Bz in axisymmetry cylindrical coordinates 
% % from input: 
% %         1. [r,z]: Position Vector investigated             [m]
% %         2. w: Angular velocity  of the rotating mirror     [rad/s]
% %         3. I: Currents of coil magnets                     [A]
% %         4. co: Coordinates of coil magnets                 [m]
% 
% %% Parameters
% g     = parameters(1);                       % Gravitational Acceleration [m/s^2]
% mu0   = parameters(2);                       % Vaccum Pemeability [N/A^2]
% rho   = parameters(3);                       % Ferrofluid material density [kg/m^3]                     
% Ms    = parameters(4);                       % Saturation Magnetization of ferrofluid [A/m]
% gamma = parameters(5); 
% 
% %% Calculating Magnetic Field H
% [Br,Bz] = multiB(r, z, I, co);       % Magnetic flux density [T]
% B       = sqrt(Br.^2+Bz.^2);         % Magnetic flux density norm [T]
% H_norm  = 1/mu0*B;                   % Magnetic field [A/m] (M = [0;0] assumed)
% 
% %% Calculate Mass Force Potentials PI
% % Langevin formula
% PI_m  = -mu0/rho*Ms*(lnsinh(H_norm * gamma)/gamma);
% PI_g  = g*z;
% PI_w  = -w^2*r.^2/2;
% PI    = PI_g+PI_w+PI_m;
% end




% function [PI, Pi_m,ratio] = massForcePotential(w,r,z,I,co)
% %massForcePotential calculates 2D mass force potential PI in axisymmetry cylindrical coordinates from input: 
% %         1. [r,z]: Position Vector investigated             [m]
% %         2. w: Angular velocity  of the rotating mirror     [rad/s]
% %         3. I: Currents of coil magnets                     [A]
% %         4. co: Coordinates of coil magnets                 [m]
% %   Some parameter values are set from this journal article:
% %   URL: https://www.researchgate.net/publication/2148158_The_Surface_Topography_of_a_Magnetic_Fluid_--_a_Quantitative_Comparison_between_Experiment_and_Numerical_Simulation
% g     = 9.806;                       % Gravitational Acceleration [m/s^2]
% mu0   = 4*pi*1e-7;                   % Vaccum Pemeability [N/A^2]
% rho   = 1236;                        % Material density [kg/m^3]
% chi0  = 1.172;                       % Magnetic susceptibility
% % chi0  = 0.01172;                     % Magnetic susceptibility
% Ms    = 14590;                       % Saturation magnetization of mateirial [A/m]          
% gamma = 3*chi0/Ms; 
% 
% [Br,Bz] = multiB(r, z, I, co);
% 
% H_norm = 1/mu0*sqrt(Br.^2+Bz.^2);    % Magnetic field [A/m] (M = [0;0] assumed)
% 
% % Pi_m  = -mu0/rho*Ms*(-(log(gamma)+log(H_norm.*csch(H_norm*gamma)))/gamma);
% Pi_m  = -mu0/rho*Ms*lnsinh(gamma * H_norm)./ gamma;
% PI    = g*z-w^2*r.^2/2+Pi_m;
% % PI = r + z;
% % PI = z;
% % ratio = Pi_m./PI
% ratio = (g*z-w^2*r.^2/2)./Pi_m;%(g*z -w^2*r.^2)/PI
% ratio = mean(abs(ratio));
% end
% 
% % PI_m  = -mu0/rho*integral(@(H_prime) M_of_Hprime(H_prime,Ms,gamma),0.00001,H_norm,'RelTol',1e-9,'AbsTol',1e-12);  % Lower bound set to .0001 rather than 0 to avoid NaN.