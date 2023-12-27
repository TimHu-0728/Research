clc
clear
close all
addpath('./used_function')

%% 2D Halbach Array Simulation
% Initialize parameters  [User Input]
w           = 1/8*0.0254;                                   % [m]     Width of magnet
h           = 1/2*0.0254;                                   % [m]     Height of magnet 
desired_len = 0.15;                                         % [m]     Desired length of the halbach array
N_coil      = ceil(desired_len/w);                          % [-]     Number of magnets
h_ff2hb     = 0.003;                                        % [m]     Ferrofluid height h
g           = 9.81;                                         % [m/s^2] Gravity 
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

rangex          = linspace(-0.1,0.1,gridsz);                             % X Range of the plot
rangey          = linspace(-0.02,0.01,gridsz);                           % Y Range of the plot
rangex_q        = linspace(-0.1,0.1,gridsz_q);                           % X Range of the plot
rangey_q        = linspace(-0.02,0.01,gridsz_q);                         % Y Range of the plot
[X,Y]           = meshgrid(rangex,rangey);
hx_eric         = (abs(rangex(end)-rangex(1)))/(gridsz - 1);             % resolution for gradient function in x dir.
hy_eric         = (abs(rangey(end)-rangey(1)))/(gridsz - 1);             % resolution for gradient function in y dir.
[Xq,Yq]         = meshgrid(rangex_q,rangey_q);

%% Calculates B fields [T]
[Bx,By]         = calculatingB(mu0,Ke,h_vec,w_vec,X,Y,theta,r0);         % B field for contour plot
[Bx_q,By_q]     = calculatingB(mu0,Ke,h_vec,w_vec,Xq,Yq,theta,r0);       % B field for quiver plot
B               = sqrt(Bx.^2+By.^2);                                     % B field in magnitude
        
%% Kelvin Body Force Map
% Ferrofluid parameters 
chi0     = 0.1;
[rho,Ms] = get_Rho_Ms(chi0);                                         %obtain rho and Ms specs from Chi0 interpolation
Ms       = 1100;                                                     % [A/m]
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
[Xq_Fk,Yq_Fk]= meshgrid(rangex_q,rangey_q);
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
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1,1);
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
levels         = linspace(0,120,20);
figure
contourf(X,Y,Fk,levels)
c              = colorbar('FontSize',22);
c.Label.String = 'Kelvin Body Force Strength (g)';
hold on
axis equal
%adding the quiver on the KBF map
quiver(Xq_Fk,Yq_Fk,Fk_x_q,Fk_y_q,'r')

% Plot the magnetized material shapes
for i = 1:size(r0,1)
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1,1);
end

title('Kelvin Body Force of Halbach Array','FontSize',22)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)

%% Plot PI tilted case
% plot the potential line for tilted and non-tilted case while maintain volume constant underneath the curve
% Theta = 0 deg (The reference case)
ff_length   = 0.05;
y_ind_above = find(rangey <= -h_ff2hb,1,'last');
X_above     = X(y_ind_above:end,:);
Y_above     = Y(y_ind_above:end,:);
Fk_above    = Fk(y_ind_above:end,:);
H_norm_above= H_norm(y_ind_above:end,:); 
rangey_above= rangey(y_ind_above:end);

delta       = 0.001;      % Ferrofluid thickness [m]
x_ff        = linspace(-ff_length,ff_length,500);
y_ff        = linspace(0,delta,200);
[X_ff,Y_ff] = meshgrid(x_ff,y_ff);

theta_0    = deg2rad(0);
PI_m       = -mu0/rho*Ms*lnsinh(gamma * H_norm_above)/ gamma;
PI_g       = g*(Y_above.* cos(theta_0) - X_above.* sin(theta_0)); 
PI         = PI_g+PI_m;
PI_0_delta = interp2(X_above,Y_above,PI,0,delta,'spline');
M          = contourc(rangex,rangey_above,PI,[PI_0_delta PI_0_delta]);
x_PI_0     = M(1,2:end);
y_PI_0     = M(2,2:end);
start_ind  = find(x_PI_0 > -ff_length,1,'first');
end_ind    = find(x_PI_0 < ff_length,1,'last');
V_ff       = trapz(x_PI_0(start_ind:end_ind),y_PI_0(start_ind:end_ind))

% Theta = 1~10 deg
theta_tilt  = deg2rad(linspace(1,10,10));

% Plot contour
figure(4)
levels2 = linspace(0,200,20);
contourf(X_above,Y_above,Fk_above,levels2)
title('PI tilted cases','FontSize',18)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)
c     = colorbar('FontSize',22);
c.Label.String = 'Kelvin Body Force Strength (g)';
hold on
patch([X_ff(1,1) X_ff(1,end) X_ff(end,end) X_ff(end,1)],[Y_ff(1,1) Y_ff(1,end) Y_ff(end,end) Y_ff(end,1)],[0.4660 0.6740 0.1880])
% Plot equipotential lines
plot(x_PI_0(start_ind:end_ind),y_PI_0(start_ind:end_ind),'r','LineWidth',1)
mkdir data
mkdir figure
writematrix([x_PI_0(start_ind:end_ind)',y_PI_0(start_ind:end_ind)'],'data/xy_0_deg.txt')
for i = 1:10
    PI_g_theta     = g*(Y_above.* cos(theta_tilt(i)) - X_above.* sin(theta_tilt(i))); 
    PI_theta       = PI_g_theta + PI_m;
    fun   = @(PI_optimal) ConservedVolume(PI_optimal,PI_theta,rangex,rangey_above,V_ff,ff_length);
    lb    = min(PI_theta,[],"all");
    ub    = max(PI_theta,[],"all");
    options = optimoptions('fmincon','Algorithm','interior-point','StepTolerance',1e-20,'OptimalityTolerance',1e-20);
    [PI_optimal,fval] = fmincon(fun,PI_0_delta,[],[],[],[],lb,ub,[],options);
    M_theta          = contourc(rangex,rangey_above,PI_theta,[PI_optimal PI_optimal]);
    x_PI_theta       = M_theta(1,2:end);
    y_PI_theta       = M_theta(2,2:end);
    start_ind_theta  = find(x_PI_theta > -ff_length,1,'first');
    end_ind_theta    = find(x_PI_theta < ff_length,1,'last');
    V_ff_theta       = trapz(x_PI_theta(start_ind_theta:end_ind_theta),y_PI_theta(start_ind_theta:end_ind_theta))
    
    % Change the equipotential line plots here!
    figure(4)
    plot(x_PI_theta(start_ind_theta:end_ind_theta),y_PI_theta(start_ind_theta:end_ind_theta),'-','LineWidth',1)
    filename = sprintf('data/xy_%d_deg.txt',i);
    writematrix([x_PI_theta(start_ind_theta:end_ind_theta)',y_PI_theta(start_ind_theta:end_ind_theta)'],filename)
    
    figure
    contourf(X_above,Y_above,PI_theta);
    filename = ['figure/',num2str(i) '_deg.png']; % Name of the file
    saveas(gcf, filename); 
end

legend('PI field(0 deg)','Ferrofluid','PI (0 deg)','PI (1 deg)','PI (2 deg)','PI (3 deg)','PI (4 deg)','PI (5 deg)','PI (6 deg)','PI (7 deg)','PI (8 deg)','PI (9 deg)','PI (10 deg)',fontsize = 16)

%% Plot Magnetization Curve from Langevin Approximation
H = linspace(0,10e6,10000);
M = Ms*(coth(gamma*H)-1/gamma./H);
figure
plot(H,M,'-',H,Ms*ones(size(H)),'--','LineWidth',1.5);
xlabel('H (A/m)','FontSize',18)
ylabel('M (A/m)','FontSize',18)
grid on
legend('M(H)','Ms',fontsize = 18)

%% Edge Effect as a function of lambda
lambda   = [1/8 1/4 1/2 1]*0.0254;                % [m]
w_ee     = lambda/4;                     % [m]
chi0     = 0.1;
[rho,Ms] = get_Rho_Ms(chi0);                                         %obtain rho and Ms specs from Chi0 interpolation
gamma    = 3*chi0/Ms;
levels   = linspace(0,120,20);
ff_length   = 0.05;
y_ind_above = find(rangey <= -h_ff2hb,1,'last');
X_above     = X(y_ind_above:end,:);
Y_above     = Y(y_ind_above:end,:);
rangey_above= rangey(y_ind_above:end);
delta       = 0.001;      % Ferrofluid thickness [m]
x_ff        = linspace(-ff_length,ff_length,500);
y_ff        = linspace(0,delta,200);
[X_ff,Y_ff] = meshgrid(x_ff,y_ff);
theta_0    = deg2rad(0);

KBF = 0;
legend_name = cell(size(w_ee));
figure
hold on
for i = 1:length(w_ee)
    N_coil = ceil(desired_len/w_ee(i));                          % [-]     Number of magnets
    if mod(N_coil,2) == 1                                       % Make even number of magnets
        N_coil = N_coil+1;
    end
    Mm          = (1.48/mu0)*ones(N_coil,1);                    % [A/m]   Magnetization
    Ke          = Mm;                                           % [A/m]   Sheet current         
    r0                          = zeros(N_coil,2);
    r0(1:N_coil/2,1)            = w_ee(i)/2;
    r0(N_coil/2+1:N_coil,1)     = -w_ee(i)/2;
    r0(1:N_coil/2,1)            = r0(1:N_coil/2,1) - w_ee(i)*(N_coil/2:-1:1)';
    r0(N_coil/2+1:N_coil,1)     = r0(N_coil/2+1:N_coil,1) + w_ee(i)*(1:1:N_coil/2)';
    r0(:,2)                     = -h/2-h_ff2hb;
    %magnet magnetization direction
    theta           = zeros(N_coil,1);
    theta(1:4:end)  = deg2rad(0);
    theta(2:4:end)  = deg2rad(90);
    theta(3:4:end)  = deg2rad(180);
    theta(4:4:end)  = deg2rad(270);
    %width input for calculating B field
    w_vec           = ones(N_coil,1)*w_ee(i);
    w_vec(2:2:end)  = h;
    
    %height input for calculating B field
    h_vec           = ones(N_coil,1)*h;
    h_vec(2:2:end)  = w_ee(i);
    
    % Calculate B field
    [Bx,By]         = calculatingB(mu0,Ke,h_vec,w_vec,X,Y,theta,r0);         % B field for contour plot
    B               = sqrt(Bx.^2+By.^2);                                     % B field in magnitude
            
    % Calculate Kelvin Body Force
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

    Fk_above    = Fk(y_ind_above:end,:);
    H_norm_above= H_norm(y_ind_above:end,:); 

    PI_m       = -mu0/rho*Ms*lnsinh(gamma * H_norm_above)/ gamma;
    PI_g       = g*(Y_above.* cos(theta_0) - X_above.* sin(theta_0)); 
    PI         = PI_g+PI_m;

    subplot(2,2,i)
    levels2 = linspace(-1,0,40);            % Set the lower bound to investigate the PI level (the white region) you want 
    contourf(X_above,Y_above,PI,levels2)
    title_name = ['PI (\lambda = ',num2str(lambda(i)/0.0254),' in.)'];
    title(title_name,'FontSize',18)
    xlabel('X (m)','FontSize',18)
    ylabel('Y (m)','FontSize',18)
    c     = colorbar('FontSize',22);
    c.Label.String = 'PI';
    hold on
    legend_name{i} = ['\lambda = ',num2str(lambda(i)/0.0254),' in.'];
end
% contour_name = ['Kelvin Body Force (',legend_name{1},')'];
% legend(contour_name,'Ferrofluid',legend_name{1},legend_name{2},legend_name{3})

%% Infinitely long Halbach Array




%% Functions
function f = ConservedVolume(PI_theta,PI,rangex,rangey,V_ff,ff_length)
    M            = contourc(rangex,rangey,PI,[PI_theta PI_theta]);
    x_PI_theta   = M(1,2:end);
    y_PI_theta   = M(2,2:end);
    start_ind    = find(x_PI_theta > -ff_length,1,'first');
    end_ind      = find(x_PI_theta < ff_length,1,'last');
    V_PI         = trapz(x_PI_theta(start_ind:end_ind),y_PI_theta(start_ind:end_ind));
    f            = abs(V_PI-V_ff);
end
