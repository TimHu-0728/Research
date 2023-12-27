% Path opreations
clc
clear
close all
addpath('./functions')

%% CONFIGURATION
% Parameters
mu    = 4*pi*1e-7;                  % Vacuum magnetic permeability (N/(A^2))
rho   = 1290;                       % Ferrofluid material density [kg/m^3]
Ms    = 28250;                      % Saturation magnetization of Ferrofluid mateirial [A/m]
chi0  = 0.1;                          % Initial Magenetic susceptibility of Ferrofluid (EMG 700 with dilution, SI unit)
gamma = 3*chi0/Ms;

% [User Input] Coils positions [x y z] (m)
Nx    = 8;
Ny    = 8;
[x,y] = meshgrid(linspace(-0.1,0.1,Nx),linspace(-0.1,0.1,Ny));
coil_positions = [y(:), x(:), zeros(Nx*Ny,1)];                  % Position of the coils

% [User Input] Coils orientation angles [psi (z-axis), theta (y-axis), phi (x-axis)] (deg), 
Euler_angles =zeros(size(coil_positions,1),3);
for i=1:size(coil_positions,1)
    Euler_angles(i,:) = deg2rad((90*(i-1)-270*(floor(i/8)))*[0 1 0]);
end
                 
% [User Input] Rectangular Coil Parameters
a  = 0.02;                     % Rectangular coil length (m)
b  = 0.02;                     % Rectangular coil width (m)
h  = 0.02;                     % Total height (m)
N  = 20;                        % Number of stacks
I  = 1.48/mu/N/h;              % Current (A)

% [User Input] Simulation Domain
xrange = [-0.11,0.11]; % Mesh range in x
yrange = [-0.11,0.11]; % in y
zrange = [1e-2,5e-2];          % in z
Nx     = 60;                  % Mesh density in x
Ny     = 60;                  % in y
Nz     = 10;                    % in z

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


% Calculating magnetic force potential PI_m for EMG 700 from Ferrotec Inc.
% (Assuming the ferrofluid's relationship between magnetization(M) and 
% magnetic field strength(H) follows a Langevin function)
% 
B                     = sqrt(Bx.^2+By.^2+Bz.^2);                      % Magnetic flux density norm [T]
H                     = 1/mu*B;                                       % Magnetic field [A/m] (M = [0;0] assumed)
PI_m                  = -mu/rho*Ms*(lnsinh(H * gamma)/gamma);         % Calculate magnetic mass force potential (Langevin curve)
[F_m_x,F_m_y,F_m_z]   = gradient(-PI_m);                              % Calculate the force vector components
level                 = 4;
F                     = sqrt(F_m_x(:,:,level).^2+F_m_y(:,:,level).^2+F_m_z(:,:,level).^2);


% Plotting Force Map
% Use quiver3 to plot B field in 3d space
scal   = 1.2;      % vector scale in quiver3
q=quiver3(x,y,z,F_m_x,F_m_y,F_m_z,scal,'LineWidth',0.75,'MaxHeadSize',1,'Color',[0.2 0.6 0.2]);
axis equal
xlabel('X (m)','FontSize',14)
ylabel('Y (m)','FontSize',14)
zlabel('Z (m)','FontSize',14)
title('Magnetic Force Map of Halbach Array (Tesla)','FontSize',16)
hold on

% plot 3d coil shape
netplot(a,b,h,N,coil_positions,Euler_angles)
hold on

sc = surfc(x(:,:,1),y(:,:,1),F*1e-6+0.02);
colorbar
sc(2).ZLocation = 'zmax';
