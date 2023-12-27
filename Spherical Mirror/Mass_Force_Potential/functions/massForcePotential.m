function PI_m = massForcePotential(w,r,z,I,co,parameters)
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
PI_m  = -mu0/rho*Ms*(-(log(gamma)+log(H_norm.*csch(H_norm*gamma)))/gamma);
% PI_g  = g*z;
% PI_w  = -w^2*r.^2/2;
% PI    = PI_g+PI_w+PI_m;
end