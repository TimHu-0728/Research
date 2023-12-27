clc
clear
close all
addpath('./functions')

%% Import Demagnetization curve data of N52 magnets
% K&J Magnetics Demagnetization curve data for N52 magents in CGS unit
% URL: https://www.kjmagnetics.com/bhcurves.asp
load("data/N52_demag_curve.mat")

%% One magnet simulation
% Perform one manet simulation to calibrate Mm so that the B field in our code matches
% the B field calculated by the K & J magnet for B888-N52 magnet
w     = 1/2*0.0254;                                   % [m]     Width
h     = 1/2*0.0254;                                   % [m]     Height
mu0   = 4*pi*1e-7;                                    % [N/A^2] Permeability of free space
N_coil= 1;   
% Calibrate Mm with K&J Magent B_field
r0    = [0 0];                                        % [m]
theta = deg2rad(0);                                   % [rad]
gridsz= 2000;                                         % Grid size for contour plot
rangex= 0.04;                                         % [m]
rangey= 0.02;                                         % [m]
B_ref = 6451.4*0.0001;                                % [T]
x_ref = 0*0.0254;                                     % [m]
y_ref = 0.25*0.0254;                                  % [m]
[Mm,B_xy] = Calibrate_Mm(B_ref,mu0,h,w,x_ref,y_ref,theta,r0);

B_ref
B_xy
Mm_vec = Mm*ones(N_coil,1);
Ke     = Mm;                                           % [A/m]   Sheet current

[X,Y] = meshgrid(linspace(-rangex,rangex,gridsz),linspace(-rangey,rangey,gridsz));
[Bx,By]     = calculatingB(mu0,Ke,h,w,X,Y,theta,r0);         % B field for contour plot
B           = sqrt(Bx.^2+By.^2);


figure
contourf(X,Y,B,linspace(0,1,15))
c  = colorbar;
axis equal
hold on
Plot_magnets(theta,w,h,r0,1)

%% 2D Halbach Array Simulation
tic
% Initialize parameters  [User Input]
% B888-N52
N_coil = 20;
mu0   = 4*pi*1e-7;                                    % [N/A^2] Permeability of free space
Mm    = (1.48/mu0)*ones(N_coil,1);                    % [A/m]   Magnetization
Ke    = Mm;                                           % [A/m]   Sheet current         
% Magnets dimensions
d     = 0.002;                                        % [m]     Gap between magnets 
w     = 1/2*0.0254;                                   % [m]     Width
h     = 1/2*0.0254;                                     % [m]     Height

% Position of the magnets
r0 = zeros(N_coil,2);
r0(1:N_coil/2,1) = w/2+d/2;
r0(N_coil/2+1:N_coil,1) = -w/2-d/2;
r0(1:N_coil/2,1) = r0(1:N_coil/2,1) - w*(N_coil/2:-1:1)'-d*(N_coil/2:-1:1)';
r0(N_coil/2+1:N_coil,1) = r0(N_coil/2+1:N_coil,1) + w*(1:1:N_coil/2)'+d*(1:1:N_coil/2)';

% Orientation of the magnets( 0: up, 90 left, 180 down, 270 right )
theta = zeros(N_coil,1);
theta(1:4:end) = deg2rad(0);
theta(2:4:end) = deg2rad(90);
theta(3:4:end) = deg2rad(180);
theta(4:4:end) = deg2rad(270);

% Width
w_vec = ones(N_coil,1)*w;
w_vec(2:2:end) = h;

% Height
h_vec = ones(N_coil,1)*h;
h_vec(2:2:end) = w;

% Sets meshgrids [User Input]
gridsz   = 2000;                                             % Grid size for contour plot
gridsz_q = 50;                                               % Grid size for quiver plot
rangex    = 1.8*max(r0(:,1));                                               % Range of the plot
rangey    = 0.8*max(r0(:,1));
[X,Y]    = meshgrid(linspace(-rangex,rangex,gridsz),linspace(-rangey,rangey,gridsz)+0.01);
hx_eric = 2*rangex/(gridsz - 1);
hy_eric = 2*rangey/(gridsz - 1);
[Xq,Yq]  = meshgrid(linspace(-rangex,rangex,gridsz_q),linspace(-rangey,rangey,gridsz_q)+0.01);

%% Calculations and plotting
% Calculates B fields [T]
[Bx,By]     = calculatingB(mu0,Ke,h_vec,w_vec,X,Y,theta,r0);         % B field for contour plot
% [Bx_q,By_q] = calculatingB(mu0,Ke,h_vec,w_vec,Xq,Yq,theta,r0);       % B field for quiver plot
B           = sqrt(Bx.^2+By.^2);



%% Find maximum B betweem the gaps (Assume coil is horizontally placed)
% ind       = zeros(gridsz,gridsz)+ (X >= (r0(1,1) - w/2 - d) & X <= (r0(1,1) - w/2) & Y <= (w/2) & Y >= (- w/2));
% for i=1:N_coil
%     ind   =  ind + (X >= (r0(i,1) + w/2) & X <= (r0(i,1) + w/2 + d) & Y <= (w/2) & Y >= (- w/2));
% end

% Find the index of space around the magnets
ind       = ones(gridsz,gridsz);
for i = 1: N_coil
    if theta(i) == 0 || theta(i) == pi
        ind = ind - (X >= (r0(i,1) - w/2) & X <= (r0(i,1) + w/2) & Y <= (r0(i,2) + h/2) & Y >= (r0(i,2) - h/2));
    elseif theta(i) == pi/2 || theta(i) == 3/2*pi
        ind = ind - (X >= (r0(i,1) - h/2) & X <= (r0(i,1) + h/2) & Y <= (r0(i,2) + w/2) & Y >= (r0(i,2) - w/2));
    end
end

B_max     = max(B(logical(ind)));                                       
ind_max   = find(abs(B-B_max)<0.000001);
X_max     = X(ind_max);
Y_max     = Y(ind_max);

% Find the maximum temperature T_max such that N52 magnets, surrounded by 
% B_max, will not be demagnetized to 96% of its original value.

T_max     = get_Tmax(-B_max,data(:,1:7),0.96);                          % [°C]

% % Plot the H_max on the 9*demagnetization curve of N52
Plotting_result(-B_max / mu0,data)

%% Kelvin Body Force Map

% Ferrofluid parameters 
chi0     = 0.1;
[rho,Ms] = get_Rho_Ms(chi0);
gamma    = 3*chi0/Ms;

Hx       = Bx / mu0;
Hy       = By / mu0;
H_norm   = sqrt(Hx.^2+Hy.^2);
M        = Ms*(coth(gamma*H_norm)-1/gamma./H_norm);
nx_H     = Hx./H_norm;
ny_H     = Hy./H_norm;
M_vec    = cat(3,M.*nx_H,M.*ny_H);

[Hx_x,Hx_y]  = gradient(Hx,hx_eric,hy_eric);
[Hy_x,Hy_y]  = gradient(Hy,hx_eric,hy_eric);

Fk_x     = mu0 * (dot(M_vec,cat(3,Hx_x,Hx_y),3));
Fk_y     = mu0 * (dot(M_vec,cat(3,Hy_x,Hy_y),3));
Fk       = sqrt(Fk_x.^2+Fk_y.^2)/rho/9.81;

% quiver mesh
gridsz_Fk     = 50;
[Xq_Fk,Yq_Fk]  =  meshgrid(linspace(-rangex,rangex,gridsz_Fk),linspace(-rangey,rangey,gridsz_Fk)+0.005);
Fk_x_q     = interp2(X,Y,Fk_x,Xq_Fk,Yq_Fk);
Fk_y_q     = interp2(X,Y,Fk_y,Xq_Fk,Yq_Fk);

%% Plot 2D Halbach Array with its B field
% Plot the B field strength contour
figure
contourf(X,Y,B,linspace(0,1,12))
c              = colorbar('FontSize',22);
c.Label.String = 'B field strength (Tesla)';
hold on
% 
% quiver(Xq,Yq,Bx_q,By_q,'Color',[1 1 1])
axis equal

% Plot the magnetized material shapes
for i = 1:size(r0,1)
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1);
end

% Plot B vector field


% Plot Ferrofluid layer
% plot(X_layer,y_layer*ones(size(X_layer)),'Color',[0.9290 0.6940 0.5000],'LineWidth',3)

grid on
title('B Field of Halbach Array','FontSize',22)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)
scatter(X_max,Y_max,40,'red','filled')
legen=legend('','','','','','','','','','','','','','','','','','','','','','','','','','','Ferrofluid layer','$B_{max}$',fontsize=16);
set(legen,'Interpreter','latex')



%% Plot the Kelvin Body Force contour map
% x_level        = linspace(0,1.01,15);
levels         = [linspace(0,120,20)];
figure
contourf(X,Y,Fk,levels)
c              = colorbar('FontSize',22);
c.Label.String = ' log_{10} (Kelvin Body Force Strength (g))';
hold on
axis equal
% quiver(Xq_Fk,Yq_Fk,Fk_x_q,Fk_y_q,'r')

% Plot the magnetized material shapes
for i = 1:size(r0,1)
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1);
end

% Plot Ferrofluid layer

title('Kelvin Body Force of Halbach Array','FontSize',22)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)



