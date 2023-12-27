
% Test For Rosensweig Instability
% An experimental study on Rosensweig instability of a ferrofluid droplet
%  Ching-Yao Chen and Z.-Y. Cheng
%  Citation: Physics of Fluids (1994-present) 20, 054105 (2008); doi: 10.1063/1.2929372
%  View online: http://dx.doi.org/10.1063/1.2929372

% Parameters
H0    = 59126.0;       % [A/m]
d_rho = 1530;          % [kg/m^3]
sigma = 25e-3;         % [N/m]
Ms    = 47746.5;       % [A/m]
chi0  = 6.79;          % 
gamma = 2*chi0/Ms;
g     = 9.806;         % m/s^2  
mu0   = 4*pi*1e-7;     % [N/A^2]

r0      = sqrt(((H0+Ms*(coth(gamma*H0)-1./(gamma*H0))).*(1+Ms*(1./(H0.^2*gamma)-gamma*csch(H0*gamma).^2)))./H0);
Mc      = sqrt(2/mu0*(1+1./r0)*sqrt(g*d_rho*sigma))              % [A/m]
M       = Ms*(coth(gamma*H0)-1/gamma./H0)                        % [A/m]





