function error = mymodel(x)
%MYMODEL This is the math model of spherical mirror used for sensitivity analysis
%--------------------------------------------------------------------------------------
%   Inputs: 
%          x1                   Density of ferrofluid
%          x2                   Saturation magnetization 
%          x3                   Initial Magenetic susceptibility
%          x4                   Angular velocity of disk
%          [x5, x6, x7, x8]     Current inside coils
%          [x9 x10;
%           x11 x12;
%           x13 x14;
%           x15 x16]            Coil positions in 2D

%---------------------------------------------------------------------------------------
%   Output: error is sum of the deviation square (at each point along the spherical 
%           mirror surface) from the ideal surface
%---------------------------------------------------------------------------------------
%   Date:   10-19-2023
%%
g                = 9.806;                   % Gravitational Acceleration [m/s^2]
mu0              = 4*pi*1e-7;               % Vaccum Pemeability [N/A^2]
rho              = x(1);                    % Ferrofluid material density [kg/m^3]
Ms               = x(2);                    % Saturation magnetization of Ferrofluid mateirial [A/m]
chi0             = x(3);                    % Initial Magenetic susceptibility of Ferrofluid (EMG 700 with dilution, SI unit)
gamma            = 3*chi0/Ms;
w                = x(4);                    % Angular velocity of disk [rad/s]
I                = [x(5);x(6);x(7);x(8)];   % Current inside each coil [A]
rz               = [x(9),x(10);    
                    x(11),x(12);
                    x(13),x(14);
                    x(15),x(16)];           % Positions of coils [m]
parameters       = [g;mu0;rho;Ms;gamma];
N_count          = 91;
r_flat           = 15;
step             = 0.0001;
[xx, yy, nx, ny] = surface_params(r_flat, N_count, step); 
X_ideal          = zeros(1, N_count);
Y_ideal          = zeros(1, N_count);
top              = massForcePotential(w,xx,yy,I,rz,parameters) ...
                   - massForcePotential(w,X_ideal,Y_ideal,I,rz,parameters);
bottom           = 1/(2*step) * (massForcePotential(w,xx+step*nx,yy+step*ny,I,rz,parameters) ...
                   - massForcePotential(w,xx-step*nx,yy-step*ny,I,rz,parameters));
error            = sum((top./bottom).^2);
end

