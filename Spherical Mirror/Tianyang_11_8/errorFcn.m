function Deviation_max = errorFcn(params)
% function y=errorFcn(g,rho,Ms,chi0,w,I1,I2,I3,I4,r1,z1,r2,z2,r3,z3,r4,z4)
%
%MYMODEL This is the math model of spherical mirror used for sensitivity analysis
%--------------------------------------------------------------------------------------
%   Inputs: 
%          g                    Gravitational acceleration       [m/s^2]
%          rho                  Density of ferrofluid            [kg/m^3]
%          Ms                   Saturation magnetization         [A/m]
%          chi0                 Initial Magenetic susceptibility
%          w                    Angular velocity of disk         [rad/s]
%          [I1, I2, I3, I4]     Current inside coils             [A]
%          [r1 z1;
%           r2 z2;
%           r3 z3;
%           r4 z4]              Coil positions in 2D             [m]

%---------------------------------------------------------------------------------------
%   Output: error is maximum of the absolute value of deviation(at each point along the spherical 
%           mirror surface) from the ideal surface
%---------------------------------------------------------------------------------------
%   Date:   10-31-2023
%%

I1 = params.DesignVars(1).Value;
I2 = params.DesignVars(2).Value;
I3 = params.DesignVars(3).Value;
I4 = params.DesignVars(4).Value;
Ms = params.DesignVars(5).Value;
chi0 = params.DesignVars(6).Value;
g = params.DesignVars(7).Value;
r1 = params.DesignVars(8).Value;
r2 = params.DesignVars(9).Value;
r3 = params.DesignVars(10).Value;
r4 = params.DesignVars(11).Value;
rho = params.DesignVars(12).Value;
w = params.DesignVars(13).Value;
z1 = params.DesignVars(14).Value;
z2 = params.DesignVars(15).Value;
z3 = params.DesignVars(16).Value;
z4 = params.DesignVars(17).Value;


mu0              = 4*pi*1e-7;               % Vaccum Pemeability [N/A^2]
gamma            = 3*chi0/Ms;
I                = [I1;I2;I3;I4];           % Current inside each coil [A]
rz               = [r1,z1;    
                    r2,z2;
                    r3,z3;
                    r4,z4];                 % Positions of coils [m]
parameters       = [g;mu0;rho;Ms;gamma];
N_count          = 91;
r_flat           = 15;
step             = 0.0001;
[xx, yy, nx, ny] = surface_params(r_flat, N_count, step); 
X_ideal          = zeros(1, N_count);
Y_ideal          = zeros(1, N_count);
top          = massForcePotential(w,xx,yy,I,rz,parameters) ...
        - massForcePotential(w,X_ideal,Y_ideal,I,rz,parameters);
bottom       = 1/(2*step) * (massForcePotential(w,xx+step*nx,yy+step*ny,I,rz,parameters) ...
        - massForcePotential(w,xx-step*nx,yy-step*ny,I,rz,parameters));
Deviation_max            = sum(abs((top./bottom))); 

% top              = massForcePotential(w,xx',yy',I,rz,parameters) ...
%                    - massForcePotential(w,X_ideal',Y_ideal',I,rz,parameters);
% bottom           = 1/(2*step) * (massForcePotential(w,(xx+step*nx)',(yy+step*ny)',I,rz,parameters) ...
%                    - massForcePotential(w,(xx-step*nx)',(yy-step*ny)',I,rz,parameters));
% Deviation_max            = sum(abs((top./bottom)));
end
