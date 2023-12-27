 function [PI,PI_m,PI_w,PI_g] = HB_massForcePotential(X_sim_range,Y_sim_range,parameters,B, th)
% massForcePotential calculates 2D mass force potentials PI, PI_m, PI_w, PI_g, 
% and magenetic flux density Br, Bz in axisymmetry cylindrical coordinates 
% from input: 
%         1. [r,z]: Position Vector investigated             [m]
%         2. w: Angular velocity  of the rotating mirror     [rad/s]
%         3. I: Currents of coil magnets                     [A]
%         4. co: Coordinates of coil magnets                 [m]
%% Parameters
g       = parameters(1);                       % Gravitational Acceleration [m/s^2]
mu0     = parameters(2);                       % Vaccum Pemeability [N/A^2]
rho     = parameters(3);                       % Ferrofluid material density [kg/m^3]                     
Ms      = parameters(4);                       % Saturation Magnetization of ferrofluid [A/m]
gamma   = parameters(5); 
sigma   = parameters(6);             % Interfacial tension                    [N/m]


% Calculating Magnetic Field H
% [B_r,B_z] = multiB(X_sim_range,Y_sim_range,mags);       % Magnetic flux density [T]
% B         = sqrt(B_r.^2+B_z.^2);                        % Magnetic flux density norm [T]
H_norm  = 1/mu0*B;                             % Magnetic field [A/m] (M = [0;0] assumed, since at the edge of surface, and M is internal only)
    
%% Calculate Mass Force Potentials PI
% Langevin formula
PI_m    = -mu0/rho*Ms*lnsinh(gamma * H_norm)./ gamma;
PI_g    = g*(Y_sim_range * cos(th) - X_sim_range * sin(th));
PI_w    = 0;
PI      = PI_g+PI_w+PI_m;
end
