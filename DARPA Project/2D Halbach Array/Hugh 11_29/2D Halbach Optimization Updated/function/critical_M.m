function [Mc,M,Hr,Hz,H0,H_norm] = critical_M(X_sim_range,Y_sim_range,parameters,n_x,n_y,mags)
%CRITICAL_M Calcultate Critical Magnetization and Magnetization of Ferrofluid at Interface.
%   This function calculates the critial magnetization Mc and magnetization M of ferrofluid at 
%   interface such that if M >= Mc, the surface of the ferrofluid will
%   have Rosensweig instability. 
%   This function takes positions, currents of the coil, and the position of
%   investigation (interface), parameters composed with Gravitational
%   acceleration g, Vacuum permeability mu0, Density difference between two
%   phase d_rho, Saturation Magnetization of ferrofluid Ms, constant gamma,
%   and Interfacial tension sigma. 

%% Extract Parameters
g       = parameters(1);             % Gravitational acceleration             [m/s^2]
mu0     = parameters(2);             % Vacuum permeability                    [N/A^2]
d_rho   = parameters(3);             % Density difference between two phase   [kg/m^3]
Ms      = parameters(4);             % Saturation Magnetization of ferrofluid [A/m]
gamma   = parameters(5);             % gamma = 3*chi0/Ms                      [m/A]
sigma   = parameters(6);             % Interfacial tension                    [N/m]

%% Calculating Normal Magnetic Field H0 Along the Interface
[B_r,B_z,B] = multiB(X_sim_range,Y_sim_range,mags);       % Magnetic flux density [T]
Hr      = B_r/mu0;
Hz      = B_z/mu0;
H0      = Hr.*n_x+Hz.*n_y;
H_norm  = 1/mu0*B;
%% Calculate r0 (geometric mean of the chord and tangent permeabilities)
r0      = sqrt(((H0+Ms*(coth(gamma*H0)-1./(gamma*H0))).*(1+Ms*(1./(H0.^2*gamma)-gamma*csch(H0*gamma).^2)))./H0);

%% Calculate Critical Magnetization Ms
Mc      = sqrt(2/mu0*(1+1./r0)*sqrt(g*d_rho*sigma));              % [A/m]
M       = Ms*(coth(gamma*H0)-1/gamma./H0);                        % [A/m]
end


