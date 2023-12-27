% Path opreations
clc
clear
close all
addpath('./functions')

%% CONFIGURATION
% Parameters
mu = 4*pi*1e-7;                % Vacuum magnetic permeability (N/(A^2))

% [User Input] Position of B field investigated (m) 
cartesian_coord=[0,0,0.1];              
x_user=cartesian_coord(1);
y_user=cartesian_coord(2);
z_user=cartesian_coord(3);

% [User Input] Coils positions [x y z] (m)
coil_positions=(-0.1:0.01:0.1)'*[1,0,0]; % Position of the coils
       

% [User Input] Coils orientation angles [psi (z-axis), theta (y-axis), phi (x-axis)] (deg), 
Euler_angles=deg2rad((0:90:21*90)'*[0,1,0]);
                 
% [User Input] Rectangular Coil Parameters
a  = 0.21;                     % Rectangular coil length (m)
b  = 0.01;                     % Rectangular coil width (m)
h  = 0.01;                     % Total height (m)
N  = 10;                       % Number of stacks
I  = 1.48/mu/N/h;              % Current (A)

% [User Input] Simulation Domain
xrange = [-0.11,0.11]; % Mesh range in x
yrange = [-0.11,0.11]; % in y
zrange = [6e-3,4e-2];  % in z
Nx     = 50; % Mesh density in x
Ny     = 50; % in y
Nz     = 20; % in z

%% COMPUTATION
% Output the magnetic field information
[Bx_user,By_user,Bz_user]=B_net(x_user,y_user,z_user,a,b,N,I,h,mu,coil_positions,Euler_angles);
fprintf('Rectangular Coil:\nLength a = %6.3f (m)\nWidth b = %6.3f (m)\nHeight h = %6.3f (m)\nNumber of Stack N = %i\nCurrent I = %6.4f (A)\n\nAt position\nr = [%6.3f %6.3f %6.3f ] (m)\n\nMagnetic Field B is\nB = [ %6.4g %6.4g %6.4g ] (Tesla)\n',a,b,h,N,I,x_user,y_user,z_user,Bx_user,By_user,Bz_user)

% Set mesh size
[x,y,z]=meshgrid(linspace(xrange(1),xrange(2),Nx),linspace(yrange(1),yrange(2),Ny),linspace(zrange(1),zrange(2),Nz));

% Calculate magnetic field at each mesh point
Bx=zeros(size(x));
By=zeros(size(y));
Bz=zeros(size(z));
for i=1:size(Bx,1)
    for j=1:size(Bx,2)
        for k=1:size(Bx,3)
            [Bx(i,j,k),By(i,j,k),Bz(i,j,k)]=B_net(x(i,j,k),y(i,j,k),z(i,j,k),a,b,N,I,h,mu,coil_positions,Euler_angles);
        end
    end
end

%% PLOTTING
% Use quiver3 to plot B field in 3d space
scal=1.2;      % vector scale in quiver3
q=quiver3(x,y,z,Bx,By,Bz,scal,'LineWidth',0.75,'MaxHeadSize',1,'Color',[0.2 0.6 0.2]);
axis equal
xlabel('X (m)','FontSize',14)
ylabel('Y (m)','FontSize',14)
zlabel('Z (m)','FontSize',14)
title('Magnetic field of Rectangular coil array (Tesla)','FontSize',16)
hold on

% plot 3d coil shape
netplot(a,b,h,N,coil_positions,Euler_angles)