clc
clear
close all

g = 9.806;
rho = 1.045639437967870e+03;
Ms = 3.054242762419179e+03;
chi0 = 0.1;
w = 0.580938597562220;
I1 = 1.276965024646607e+07;
I2 = 4.635658998444038e+06;
I3 = 2.598128595057578;
I4 = 3.480704452324897;
r1 = 18.023533723368903;
z1 = 3.213780836995019;
r2 = 7.170985872272569;
z2 = -5.999917115991651;
r3 = 24.999999968462504;
z3 = -5.999999857702986;
r4 = 24.999999998194050;
z4 = -5.999999794564147;



error = mymodel(g,rho,Ms,chi0,w,I1,I2,I3,I4,r1,z1,r2,z2,r3,z3,r4,z4)


function error = mymodel(g,rho,Ms,chi0,w,I1,I2,I3,I4,r1,z1,r2,z2,r3,z3,r4,z4)
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
top              = massForcePotential(w,xx',yy',I,rz,parameters) ...
                   - massForcePotential(w,X_ideal',Y_ideal',I,rz,parameters);
bottom           = 1/(2*step) * (massForcePotential(w,(xx+step*nx)',(yy+step*ny)',I,rz,parameters) ...
                   - massForcePotential(w,(xx-step*nx)',(yy-step*ny)',I,rz,parameters));
error            = max(abs((top./bottom)));
end