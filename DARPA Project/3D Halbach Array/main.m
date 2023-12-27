clc;clear;close all
addpath('functions/','data/')
format default

%%

tic
mu0        = 4*pi*1e-7;
% N52 B666 magnets
Mm         = 1.48/mu0;                        % [A/m]
a          = 1/8*0.0254;                         % [m]
b          = 1/8*0.0254;
h          = 1*0.0254;
dimensions = [a, b, h];                       % [m]

xi         = -0.04;
xf         = 0.04;
yi         = -0.04;
yf         = 0.04;
zi         = -0.02;
zf         = 0.02;

% Simulation field range [m]
[x,y,z]    = Mesh(2*xi,2*xf,2*yi,2*yf,2*zi,2*zf,75);          % [m]
[X,Y,Z]    = Mesh(xi,xf,yi,yf,zi,zf,30);

% One magnet
% r_mag      = [0 0 0];
% orient_mag = deg2rad([0 0 0]);
% 
% [Bx,By,Bz] = B_multiple(x,y,z,X,Y,Z,Mm,dimensions,r_mag,orient_mag,"spline");

% Halbach Array
HB_range      = 0.03;                         % [m]
n_mag      = 5;
[r_mag,orient_mag] = magnetPose(n_mag,HB_range,-0.002);

[Bx,By,Bz] = B_multiple(x,y,z,X,Y,Z,Mm,dimensions,r_mag,orient_mag,"spline");
%%
plot_HBArray(X,Y,Z,Bx,By,Bz,r_mag,orient_mag,dimensions,1)

toc



