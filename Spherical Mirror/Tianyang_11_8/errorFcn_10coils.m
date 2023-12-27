function Deviation_max = errorFcn_10coils(x)
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

w=x(1);
I1=x(2);
z1=x(3);
r1=x(4);
I2=x(5);
z2=x(6);
r2=x(7);
I3=x(8);
z3=x(9);
r3=x(10);
I4=x(11);
z4=x(12);
r4=x(13);
I5=x(14);
z5=x(15);
r5=x(16);
I6=x(17);
z6=x(18);
r6=x(19);
I7=x(20);
z7=x(21);
r7=x(22);
I8=x(23);
z8=x(24);
r8=x(25);
I9=x(26);
z9=x(27);
r9=x(28);
I10=x(29);
z10=x(30);
r10=x(31);
g=x(32);
mu0=x(33);
rho=x(34);
Ms=x(35);
chi0=x(36);
sigma=x(37);

gamma= 3*chi0/Ms;
I                = [I1;I2;I3;I4;I5;I6;I7;I8;I9;I10];           % Current inside each coil [A]
rz               = [r1,z1;    
                    r2,z2;
                    r3,z3;
                    r4,z4;
                    r5,z5;
                    r6,z6;
                    r7,z7;
                    r8,z8;
                    r9,z9;
                    r10,z10];                 % Positions of coils [m]
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
Deviation_max            = max(abs((top./bottom))); 

% top              = massForcePotential(w,xx',yy',I,rz,parameters) ...
%                    - massForcePotential(w,X_ideal',Y_ideal',I,rz,parameters);
% bottom           = 1/(2*step) * (massForcePotential(w,(xx+step*nx)',(yy+step*ny)',I,rz,parameters) ...
%                    - massForcePotential(w,(xx-step*nx)',(yy-step*ny)',I,rz,parameters));
% Deviation_max            = sum(abs((top./bottom)));
end
