clc
clear
close all
addpath('./used_function')

%% 2D Halbach Array Simulation
% Initialize parameters  [User Input]
w           = 1/8*0.0254;                                   % [m]     Width of magnet
h           = 1*0.0254;                                     % [m]     Height of magnet 
desired_len = 0.15;                                         % [m]     Desired length of the halbach array
N_coil      = ceil(desired_len/w);                                           % number of magnets
h_ff2hb     = 0.008;
g           = 9.81;                                          % [m/s^2]
mu0         = 4*pi*1e-7;                                    % [N/A^2] Permeability of free space
Mm          = (1.48/mu0)*ones(N_coil,1);                    % [A/m]   Magnetization
Ke          = Mm;                                           % [A/m]   Sheet current         

% Defining simulation resolution [User Input]
gridsz      = 3000;                                         % Grid size for contour plot
gridsz_q    = 50;                                           % Grid size for quiver plot

%magnet placement (x,y coordinates)
r0                          = zeros(N_coil,2);
r0(1:N_coil/2,1)            = w/2;
r0(N_coil/2+1:N_coil,1)     = -w/2;
r0(1:N_coil/2,1)            = r0(1:N_coil/2,1) - w*(N_coil/2:-1:1)';
r0(N_coil/2+1:N_coil,1)     = r0(N_coil/2+1:N_coil,1) + w*(1:1:N_coil/2)';
r0(:,2)                     = -h/2-h_ff2hb;
%magnet magnetization direction
theta           = zeros(N_coil,1);
theta(1:4:end)  = deg2rad(0);
theta(2:4:end)  = deg2rad(90);
theta(3:4:end)  = deg2rad(180);
theta(4:4:end)  = deg2rad(270);

%width input for calculating B field
w_vec           = ones(N_coil,1)*w;
w_vec(2:2:end)  = h;

%height input for calculating B field
h_vec           = ones(N_coil,1)*h;
h_vec(2:2:end)  = w;

rangex          = linspace(-0.1,0.1,gridsz);                                        % X Range of the plot
rangey          = linspace(-0.04,0.02,gridsz);                                      % Y Range of the plot
rangex_PI       = linspace(-0.1,0.1,gridsz);                                        % X Range of the plot
rangey_PI       = linspace(-0.04,0.02,gridsz);                                      % Y Range of the plot
rangex_q        = linspace(-0.1,0.1,gridsz_q);                                      % X Range of the plot
rangey_q        = linspace(-0.04,0.01,gridsz_q);                                    % Y Range of the plot
[X,Y]           = meshgrid(rangex,rangey);
hx_eric         = (abs(rangex(end)-rangex(1)))/(gridsz - 1);                                 % resolution for gradient function in x dir.
hy_eric         = (abs(rangey(end)-rangey(1)))/(gridsz - 1);                                 % resolution for gradient function in y dir.
[Xq,Yq]         = meshgrid(rangex_q,rangey_q);

%% Calculates B fields [T]
[Bx,By]         = calculatingB(mu0,Ke,h_vec,w_vec,X,Y,theta,r0);         % B field for contour plot
[Bx_q,By_q]     = calculatingB(mu0,Ke,h_vec,w_vec,Xq,Yq,theta,r0);       % B field for quiver plot
B               = sqrt(Bx.^2+By.^2);                                     % B field in magnitude
        
%% Kelvin Body Force Map
% Ferrofluid parameters 
chi0     = 0.1;
[rho,Ms] = get_Rho_Ms(chi0);                                         %obtain rho and Ms specs from Chi0 interpolation
gamma    = 3*chi0/Ms;

Hx       = Bx / mu0;                                                 %obtaining X magnetic field from B field
Hy       = By / mu0;                                                 %obtaining Y magnetic field from B field
H_norm   = sqrt(Hx.^2+Hy.^2);
M        = Ms*(coth(gamma*H_norm)-1/gamma./H_norm);                  %M magnitude
nx_H     = Hx./H_norm;                                               %unit vector of H in x dir
ny_H     = Hy./H_norm;                                               %unit vector of H in y dir
M_vec    = cat(3,M.*nx_H,M.*ny_H);                                   %M as vector form, as M is in direction of H

[Hx_x,Hx_y]  = gradient(Hx,hx_eric,hy_eric);                         %intermediate step for calculation of gradient
[Hy_x,Hy_y]  = gradient(Hy,hx_eric,hy_eric);                         %intermediate step for calculation of gradient   

Fk_x         = mu0 * (dot(M_vec,cat(3,Hx_x,Hx_y),3));                %Kelvin body force in x dir
Fk_y         = mu0 * (dot(M_vec,cat(3,Hy_x,Hy_y),3));                %Kelvin body force in y dir  
Fk           = sqrt(Fk_x.^2+Fk_y.^2)/rho/9.81;                       %[m/s^2] Kelvin body force magnitude, in g

% quiver mesh
[Xq_Fk,Yq_Fk]=  meshgrid(rangex_q,rangey_q);
Fk_x_q       = interp2(X,Y,Fk_x,Xq_Fk,Yq_Fk);                          %interp between known kelvin forces in x dir
Fk_y_q       = interp2(X,Y,Fk_y,Xq_Fk,Yq_Fk);                          %interp between known kelvin forces in y dir

%% Plot 2D Halbach Array with its B field
% Plot the B field strength contour
figure
contourf(X,Y,B,linspace(0,1,12))
c              = colorbar('FontSize',22);
c.Label.String = 'B field strength (Tesla)';
hold on

quiver(Xq,Yq,Bx_q,By_q,'Color',[1 1 1])
axis equal

% Plot the magnetized material shapes
for i = 1:size(r0,1)
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1);
end

% Plot Ferrofluid layer
% plot(X_layer,y_layer*ones(size(X_layer)),'Color',[0.9290 0.6940 0.5000],'LineWidth',3)

% Plot B vector field
grid on
title('B Field of Halbach Array','FontSize',22)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)
% scatter(X_max,Y_max,40,'red','filled')
% legen=legend('','','','','','','','','','','','','','','Ferrofluid layer','$B_{max}$',fontsize=16);
% set(legen,'Interpreter','latex')

%% Plot the Kelvin Body Force contour map
% x_level        = linspace(0,1.01,15);
levels         = [linspace(0,120,20)];
figure
contourf(X,Y,Fk,levels)
c              = colorbar('FontSize',22);
c.Label.String = 'Kelvin Body Force Strength (g)';
hold on

%adding the quiver on the KBF map
quiver(Xq_Fk,Yq_Fk,Fk_x_q,Fk_y_q,'r')

% Plot the magnetized material shapes
for i = 1:size(r0,1)
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1);
end

title('Kelvin Body Force of Halbach Array','FontSize',22)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)

%% Plot PI tilted case
%plot the potential line for tilted and non-tilted case while maintain volume constant underneath the curve




h_ff2hb     = 0.008;
delta       = 0.002;      % Ferrofluid thickness [m]
xrange      = max(r0(:,1))+w/2;
x_ff        = linspace(-xrange,xrange,500);
y_ff        = linspace(0,delta,200);
[X_ff,Y_ff] = meshgrid(x_ff,y_ff);

% Theta = 0 deg (The reference case)
theta_0  = deg2rad(0);
PI_m     = -mu0/rho*Ms*lnsinh(gamma * H_norm)/ gamma;
PI_g     = g*(Y.* cos(theta_0) - X.* sin(theta_0)); 
PI       = PI_g+PI_m;
PI_0_delta = interp2(X,Y,PI,0,delta,'spline');
M        = contourc(rangex,rangey,PI,[PI_0_delta PI_0_delta]);
x_PI_0   = M(1,2:end);
y_PI_0   = M(2,2:end);
start_ind = find(abs(x_PI_0 - (min(r0(:,1)) - w/2)) <= 0.0001,1,'first');
end_ind   = find(abs(x_PI_0 - (max(r0(:,1)) + w/2)) <= 0.0001,1,'last');
V_ff     = trapz(x_PI_0(start_ind:end_ind),y_PI_0(start_ind:end_ind))

% Theta = 1~10 deg
theta_tilt  = deg2rad(linspace(1,10,10));

% Plot contour
figure
levels = linspace(-5,1,10);
contourf(X,Y,PI,levels)
title('PI tilted cases','FontSize',18)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)
colorbar
hold on
patch([X_ff(1,1) X_ff(1,end)  X_ff(end,end) X_ff(end,1)],[Y_ff(1,1) Y_ff(1,end) Y_ff(end,end) Y_ff(end,1)],[0.4660 0.6740 0.1880])
% Plot equipotential lines
plot(x_PI_0(start_ind:end_ind),y_PI_0(start_ind:end_ind),'r','LineWidth',1)
mkdir xydata
writematrix([x_PI_0(start_ind:end_ind)',y_PI_0(start_ind:end_ind)'],'xydata/xy_0_deg.txt')
for i = 1:10
    PI_g_theta     = g*(Y.* cos(theta_tilt(i)) - X.* sin(theta_tilt(i))); 
    PI_theta       = PI_g_theta + PI_m;
    fun   = @(PI_optimal) ConservedVolume(PI_optimal,PI_theta,rangex,rangey,r0,w,V_ff);
    lb    = min(PI_theta,[],"all");
    ub    = max(PI_theta,[],"all");
    options = optimoptions('fmincon','Algorithm','interior-point','StepTolerance',1e-20,'OptimalityTolerance',1e-20);
    [PI_optimal,fval] = fmincon(fun,PI_0_delta,[],[],[],[],lb,ub,[],options);
    M_theta          = contourc(rangex,rangey,PI_theta,[PI_optimal PI_optimal]);
    x_PI_theta       = M_theta(1,2:end);
    y_PI_theta       = M_theta(2,2:end);
    start_ind_theta  = find(abs(x_PI_theta - (min(r0(:,1)) - w/2)) <= 0.00001,1,'first');
    end_ind_theta    = find(abs(x_PI_theta - (max(r0(:,1)) + w/2)) <= 0.00001,1,'last');
    V_ff_theta       = trapz(x_PI_theta(start_ind_theta:end_ind_theta),y_PI_theta(start_ind_theta:end_ind_theta))
    
    % Change the equipotential line plots here!
    plot(x_PI_theta(start_ind_theta:end_ind_theta),y_PI_theta(start_ind_theta:end_ind_theta),'-','LineWidth',1)
    filename = sprintf('xydata/xy_%d_deg.txt',i);
    writematrix([x_PI_theta(start_ind_theta:end_ind_theta)',y_PI_theta(start_ind_theta:end_ind_theta)'],filename)
end

% Plot the magnetized material shapes
for i = 1:size(r0,1)
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1);
end

legend('PI field(0 deg)','Ferrofluid','PI (0 deg)','PI (1 deg)','PI (2 deg)','PI (3 deg)','PI (4 deg)','PI (5 deg)','PI (6 deg)','PI (7 deg)','PI (8 deg)','PI (9 deg)','PI (10 deg)',fontsize = 16)

%% Functions
function f = ConservedVolume(PI_theta,PI,rangex,rangey,r0,w,V_ff)
    M            = contourc(rangex,rangey,PI,[PI_theta PI_theta]);
    x_PI_theta   = M(1,2:end);
    y_PI_theta   = M(2,2:end);
    start_ind    = find(abs(x_PI_theta - (min(r0(:,1)) - w/2)) <= 0.0001,1,'first');
    end_ind      = find(abs(x_PI_theta - (max(r0(:,1)) + w/2)) <= 0.0001,1,'last');
    V_PI         = trapz(x_PI_theta(start_ind:end_ind),y_PI_theta(start_ind:end_ind));
    f            = abs(V_PI-V_ff);
end
