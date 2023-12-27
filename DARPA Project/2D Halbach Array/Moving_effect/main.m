clc
clear
close all
addpath('./functions','./figure')

%%
% Initialize parameters  [User Input]
tic

lambda      = 2*0.0254;
w           = lambda/4;                                     % [m]     Width of magnet
h           = lambda/4;                                     % [m]     Height of magnet
d           = w/4;                                          % [m]     Distance between magnets
desired_len = [0.5 0.25 0.15];                              % [m]     Desired length of the halbach array
gridsz      = 1000;                                         % Grid size for contour plot, increase this to get more accuracy but need more computation cost
h_ff2hb     = 0.004;                                        % [m]     Ferrofluid height h
g           = 9.81;                                         % [m/s^2] Gravity
mu0         = 4*pi*1e-7;
N_coil      = zeros(3,3);
chi0        = 0.181;
[rho,Ms,phi]= get_Rho_Ms_Phi(chi0);                                         %obtain rho and Ms specs from Chi0 interpolation
gamma       = 3*chi0/Ms;

for i = 1:length(desired_len)
    for j = 1:length(w)
        N_coil(i,j)   = ceil(desired_len(i)/w(j));                          % [-]     Number of magnets
        if mod(N_coil(i,j),2) == 1                                       % Make even number of magnets
            N_coil(i,j) = N_coil(i,j)+1;
        end
    end
end
rangex          = linspace(-0.3,0.3,gridsz);                             % X Range of the plot
rangey          = linspace(-0.1,0.1,gridsz);                           % Y Range of the plot
[X,Y]           = meshgrid(rangex,rangey);
hx_eric         = (abs(rangex(end)-rangex(1)))/(gridsz - 1);             % resolution for gradient function in x dir.
hy_eric         = (abs(rangey(end)-rangey(1)))/(gridsz - 1);             % resolution for gradient function in y dir.


%% Calculate infinite halbach array PI
N_coil_inf           = 500;
Mm                   = (1.48/mu0)*ones(N_coil_inf,1);                    % [A/m]   Magnetization
Mm(N_coil_inf/2+1,1) = 1*Mm(N_coil_inf/2+1,1);                        % [A/m]   Change the magnetization of one magnet by 1%
Ke                   = Mm;                                               % [A/m]   Sheet current



%magnet placement (x,y coordinates)
r0                              = zeros(N_coil_inf,2);
r0(1:N_coil_inf/2,1)            = w/2+d/2;
r0(N_coil_inf/2+1:N_coil_inf,1) = -w/2-d/2;
r0(1:N_coil_inf/2,1)            = r0(1:N_coil_inf/2,1) - w*(N_coil_inf/2:-1:1)' - d*(N_coil_inf/2:-1:1)';
r0(N_coil_inf/2+1:N_coil_inf,1) = r0(N_coil_inf/2+1:N_coil_inf,1) + w*(1:1:N_coil_inf/2)' + d*(1:1:N_coil_inf/2)';
r0(:,2)                         = -h/2-h_ff2hb;

%magnet magnetization direction
theta0           = zeros(N_coil_inf,1);
theta0(1:4:end)  = deg2rad(0);
theta0(2:4:end)  = deg2rad(90);
theta0(3:4:end)  = deg2rad(180);
theta0(4:4:end)  = deg2rad(270);

%width input for calculating B field
w_vec0           = ones(N_coil_inf,1)*w;
w_vec0(2:2:end)  = h;

%height input for calculating B field
h_vec0           = ones(N_coil_inf,1)*h;
h_vec0(2:2:end)  = w;

y_ind_above = 300;
y_ind_below = 700;
% Calculates B fields [T]
[Bx,By]         = calculatingB(mu0,Ke,h_vec0,w_vec0,X,Y,theta0,r0);         % B field for contour plot
B               = sqrt(Bx.^2+By.^2);
Hx              = Bx / mu0;                                                 %obtaining X magnetic field from B field
Hy              = By / mu0;                                                 %obtaining Y magnetic field from B field
H_norm          = sqrt(Hx.^2+Hy.^2);


% PI_ref
X_above     = X(y_ind_above:y_ind_below,:);
Y_above     = Y(y_ind_above:y_ind_below,:);
H_norm_above= H_norm(y_ind_above:y_ind_below,:);
theta_0     = deg2rad(0);
PI_m        = -mu0/rho*Ms*lnsinh(gamma * H_norm_above)/ gamma;
PI_g        = g*(Y_above.* cos(theta_0) - X_above.* sin(theta_0));
PI_ref      = PI_g+PI_m;
% Gradient of PI_ref
[nabla_PI_ref_x,nabla_PI_ref_y] = gradient(PI_ref, hx_eric, hy_eric);

% Plot PI_inf Contour
% figure
% %subplot(2,2,j)
% levels1 = linspace(-2,0,15);
% contourf(X_above,Y_above,PI_ref,levels1)
% titlename = 'PI_{ref}';
% title(titlename,'FontSize',18)
% xlabel('X (m)','FontSize',18)
% ylabel('Y (m)','FontSize',18)
% c         = colorbar('FontSize',18);
% ylabel(c,'PI [m^2/s^2]','FontSize',18);
% hold on
% % % Plot the magnetized material shapes
% for k = 1:size(r0,1)
%     Plot_magnets(theta0(k),w_vec0(k),h_vec0(k),r0(k,:),1,1,N_coil_inf,N_coil_inf);
% end
% xlim([-0.3 0.3])
% axis equal

%% Calculate infinite halbach array with 1 magnet moving PI_moved 1000 um
delta  = 1000e-6;
y_move = [delta -delta 0 0]';
x_move = [0 0 delta -delta]';

for i = 1:4

%magnet placement (x,y coordinates)
r0                              = zeros(N_coil_inf,2);
r0(1:N_coil_inf/2,1)            = w/2+d/2;
r0(N_coil_inf/2+1:N_coil_inf,1) = -w/2-d/2;
r0(1:N_coil_inf/2,1)            = r0(1:N_coil_inf/2,1) - w*(N_coil_inf/2:-1:1)' - d*(N_coil_inf/2:-1:1)';
r0(N_coil_inf/2+1:N_coil_inf,1) = r0(N_coil_inf/2+1:N_coil_inf,1) + w*(1:1:N_coil_inf/2)' + d*(1:1:N_coil_inf/2)';
r0(:,2)                         = -h/2-h_ff2hb;
    
r0(N_coil_inf/2+1,2) = r0(N_coil_inf/2+1,2)+y_move(i);
r0(N_coil_inf/2+1,1) = r0(N_coil_inf/2+1,1)+x_move(i);

%magnet magnetization direction
theta0           = zeros(N_coil_inf,1);
theta0(1:4:end)  = deg2rad(0);
theta0(2:4:end)  = deg2rad(90);
theta0(3:4:end)  = deg2rad(180);
theta0(4:4:end)  = deg2rad(270);

%width input for calculating B field
w_vec0           = ones(N_coil_inf,1)*w;
w_vec0(2:2:end)  = h;

%height input for calculating B field
h_vec0           = ones(N_coil_inf,1)*h;
h_vec0(2:2:end)  = w;

y_ind_above = 300;
y_ind_below = 700;
% Calculates B fields [T]
[Bx,By]         = calculatingB(mu0,Ke,h_vec0,w_vec0,X,Y,theta0,r0);         % B field for contour plot
B               = sqrt(Bx.^2+By.^2);
Hx              = Bx / mu0;                                                 %obtaining X magnetic field from B field
Hy              = By / mu0;                                                 %obtaining Y magnetic field from B field
H_norm          = sqrt(Hx.^2+Hy.^2);


% PI_ref
X_above     = X(y_ind_above:y_ind_below,:);
Y_above     = Y(y_ind_above:y_ind_below,:);
H_norm_above= H_norm(y_ind_above:y_ind_below,:);
theta_0     = deg2rad(0);
PI_m        = -mu0/rho*Ms*lnsinh(gamma * H_norm_above)/ gamma;
PI_g        = g*(Y_above.* cos(theta_0) - X_above.* sin(theta_0));
PI_moved    = PI_g+PI_m;

% Calculating moving effect
eps = abs((PI_moved-PI_ref)./nabla_PI_ref_y)*1e6;

% Plot PI_moved
% figure
% %subplot(2,2,j)
% levels1 = linspace(-2,0,15);
% contourf(X_above,Y_above,PI_moved,levels1)
% titlename = ['PI_{moved} ({\delta}y = ', num2str(y_move(i)*1e6),' {\mu}m, {\delta}x = ', num2str(x_move(i)*1e6),' {\mu}m)'];
% title(titlename,'FontSize',18)
% xlabel('X (m)','FontSize',18)
% ylabel('Y (m)','FontSize',18)
% c         = colorbar('FontSize',18);
% ylabel(c,'PI [m^2/s^2]','FontSize',18);
% hold on
% % % Plot the magnetized material shapes
% for k = 1:size(r0,1)
%     Plot_magnets(theta0(k),w_vec0(k),h_vec0(k),r0(k,:),1,1,k,N_coil_inf);
% end
% 
% xlim([-0.3 0.3])
% axis equal

% Plot epsilon
figure
%subplot(2,2,j)
levels1 = linspace(-0.002,0.002,20);
contourf(X_above,Y_above,log10(eps),[-5,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.5,4,4.5,5,5.5,6])
titlename = ['|\epsilon| ({\delta}y = ', num2str(y_move(i)*1e6),' {\mu}m, {\delta}x = ', num2str(x_move(i)*1e6),' {\mu}m)'];
title(titlename,'FontSize',24)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)
c         = colorbar('FontSize',18);
ylabel(c,'log_{10}(|\epsilon| [{\mu}m])','FontSize',24);
hold on
% % Plot the magnetized material shapes
for k = 1:size(r0,1)
    Plot_magnets(theta0(k),w_vec0(k),h_vec0(k),r0(k,:),1,1,k,N_coil_inf);
end
xlim([-0.1 0.1])
axis equal

if i == 1
figname = 'figure/y + 1000um.fig';
elseif i == 2
figname = 'figure/y - 1000um.fig';
elseif i == 3
figname = 'figure/x + 1000um.fig';
elseif i == 4
figname = 'figure/x - 1000um.fig';
end
saveas(gcf,figname)
end

%% Calculate infinite halbach array with 1 magnet moving PI_moved 15 um
delta  = 15e-6;
y_move = [delta -delta 0 0]';
x_move = [0 0 delta -delta]';

for i = 1:4

%magnet placement (x,y coordinates)
r0                              = zeros(N_coil_inf,2);
r0(1:N_coil_inf/2,1)            = w/2+d/2;
r0(N_coil_inf/2+1:N_coil_inf,1) = -w/2-d/2;
r0(1:N_coil_inf/2,1)            = r0(1:N_coil_inf/2,1) - w*(N_coil_inf/2:-1:1)' - d*(N_coil_inf/2:-1:1)';
r0(N_coil_inf/2+1:N_coil_inf,1) = r0(N_coil_inf/2+1:N_coil_inf,1) + w*(1:1:N_coil_inf/2)' + d*(1:1:N_coil_inf/2)';
r0(:,2)                         = -h/2-h_ff2hb;
    
r0(N_coil_inf/2+1,2) = r0(N_coil_inf/2+1,2)+y_move(i);
r0(N_coil_inf/2+1,1) = r0(N_coil_inf/2+1,1)+x_move(i);

%magnet magnetization direction
theta0           = zeros(N_coil_inf,1);
theta0(1:4:end)  = deg2rad(0);
theta0(2:4:end)  = deg2rad(90);
theta0(3:4:end)  = deg2rad(180);
theta0(4:4:end)  = deg2rad(270);

%width input for calculating B field
w_vec0           = ones(N_coil_inf,1)*w;
w_vec0(2:2:end)  = h;

%height input for calculating B field
h_vec0           = ones(N_coil_inf,1)*h;
h_vec0(2:2:end)  = w;

y_ind_above = 300;
y_ind_below = 700;
% Calculates B fields [T]
[Bx,By]         = calculatingB(mu0,Ke,h_vec0,w_vec0,X,Y,theta0,r0);         % B field for contour plot
B               = sqrt(Bx.^2+By.^2);
Hx              = Bx / mu0;                                                 %obtaining X magnetic field from B field
Hy              = By / mu0;                                                 %obtaining Y magnetic field from B field
H_norm          = sqrt(Hx.^2+Hy.^2);


% PI_ref
X_above     = X(y_ind_above:y_ind_below,:);
Y_above     = Y(y_ind_above:y_ind_below,:);
H_norm_above= H_norm(y_ind_above:y_ind_below,:);
theta_0     = deg2rad(0);
PI_m        = -mu0/rho*Ms*lnsinh(gamma * H_norm_above)/ gamma;
PI_g        = g*(Y_above.* cos(theta_0) - X_above.* sin(theta_0));
PI_moved    = PI_g+PI_m;

% Calculating moving effect
eps = abs((PI_moved-PI_ref)./nabla_PI_ref_y)*1e6;

% Plot PI_moved
% figure
% %subplot(2,2,j)
% levels1 = linspace(-2,0,15);
% contourf(X_above,Y_above,PI_moved,levels1)
% titlename = ['PI_{moved} ({\delta}y = ', num2str(y_move(i)*1e6),' {\mu}m, {\delta}x = ', num2str(x_move(i)*1e6),' {\mu}m)'];
% title(titlename,'FontSize',18)
% xlabel('X (m)','FontSize',18)
% ylabel('Y (m)','FontSize',18)
% c         = colorbar('FontSize',18);
% ylabel(c,'PI [m^2/s^2]','FontSize',18);
% hold on
% % % Plot the magnetized material shapes
% for k = 1:size(r0,1)
%     Plot_magnets(theta0(k),w_vec0(k),h_vec0(k),r0(k,:),1,1,k,N_coil_inf);
% end
% 
% xlim([-0.3 0.3])
% axis equal

% Plot epsilon
figure
%subplot(2,2,j)
levels1 = linspace(-0.00002,0.00002,20);
contourf(X_above,Y_above,log10(eps),[-6,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4])
titlename = ['|\epsilon| ({\delta}y = ', num2str(y_move(i)*1e6),' {\mu}m, {\delta}x = ', num2str(x_move(i)*1e6),' {\mu}m)'];
title(titlename,'FontSize',24)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)
c         = colorbar('FontSize',18);
ylabel(c,'log_{10}(|\epsilon| [{\mu}m])','FontSize',24);
hold on
% % Plot the magnetized material shapes
for k = 1:size(r0,1)
    Plot_magnets(theta0(k),w_vec0(k),h_vec0(k),r0(k,:),1,1,k,N_coil_inf);
end
xlim([-0.1 0.1])
axis equal

if i == 1
figname = 'figure/y + 15um.fig';
elseif i == 2
figname = 'figure/y - 15um.fig';
elseif i == 3
figname = 'figure/x + 15um.fig';
elseif i == 4
figname = 'figure/x - 15um.fig';
end
saveas(gcf,figname)
end

%% Calculate infinite halbach array with 1 magnet moving PI_moved 1 um
delta  = 1e-6;
y_move = [delta -delta 0 0]';
x_move = [0 0 delta -delta]';

for i = 1:4

%magnet placement (x,y coordinates)
r0                              = zeros(N_coil_inf,2);
r0(1:N_coil_inf/2,1)            = w/2+d/2;
r0(N_coil_inf/2+1:N_coil_inf,1) = -w/2-d/2;
r0(1:N_coil_inf/2,1)            = r0(1:N_coil_inf/2,1) - w*(N_coil_inf/2:-1:1)' - d*(N_coil_inf/2:-1:1)';
r0(N_coil_inf/2+1:N_coil_inf,1) = r0(N_coil_inf/2+1:N_coil_inf,1) + w*(1:1:N_coil_inf/2)' + d*(1:1:N_coil_inf/2)';
r0(:,2)                         = -h/2-h_ff2hb;
    
r0(N_coil_inf/2+1,2) = r0(N_coil_inf/2+1,2)+y_move(i);
r0(N_coil_inf/2+1,1) = r0(N_coil_inf/2+1,1)+x_move(i);

%magnet magnetization direction
theta0           = zeros(N_coil_inf,1);
theta0(1:4:end)  = deg2rad(0);
theta0(2:4:end)  = deg2rad(90);
theta0(3:4:end)  = deg2rad(180);
theta0(4:4:end)  = deg2rad(270);

%width input for calculating B field
w_vec0           = ones(N_coil_inf,1)*w;
w_vec0(2:2:end)  = h;

%height input for calculating B field
h_vec0           = ones(N_coil_inf,1)*h;
h_vec0(2:2:end)  = w;

y_ind_above = 300;
y_ind_below = 700;
% Calculates B fields [T]
[Bx,By]         = calculatingB(mu0,Ke,h_vec0,w_vec0,X,Y,theta0,r0);         % B field for contour plot
B               = sqrt(Bx.^2+By.^2);
Hx              = Bx / mu0;                                                 %obtaining X magnetic field from B field
Hy              = By / mu0;                                                 %obtaining Y magnetic field from B field
H_norm          = sqrt(Hx.^2+Hy.^2);


% PI_ref
X_above     = X(y_ind_above:y_ind_below,:);
Y_above     = Y(y_ind_above:y_ind_below,:);
H_norm_above= H_norm(y_ind_above:y_ind_below,:);
theta_0     = deg2rad(0);
PI_m        = -mu0/rho*Ms*lnsinh(gamma * H_norm_above)/ gamma;
PI_g        = g*(Y_above.* cos(theta_0) - X_above.* sin(theta_0));
PI_moved    = PI_g+PI_m;

% Calculating moving effect
eps = abs((PI_moved-PI_ref)./nabla_PI_ref_y)*1e6;

% Plot PI_moved
% figure
% %subplot(2,2,j)
% levels1 = linspace(-2,0,15);
% contourf(X_above,Y_above,PI_moved,levels1)
% titlename = ['PI_{moved} ({\delta}y = ', num2str(y_move(i)*1e6),' {\mu}m, {\delta}x = ', num2str(x_move(i)*1e6),' {\mu}m)'];
% title(titlename,'FontSize',18)
% xlabel('X (m)','FontSize',18)
% ylabel('Y (m)','FontSize',18)
% c         = colorbar('FontSize',18);
% ylabel(c,'PI [m^2/s^2]','FontSize',18);
% hold on
% % % Plot the magnetized material shapes
% for k = 1:size(r0,1)
%     Plot_magnets(theta0(k),w_vec0(k),h_vec0(k),r0(k,:),1,1,k,N_coil_inf);
% end
% 
% xlim([-0.3 0.3])
% axis equal

% Plot epsilon
figure
%subplot(2,2,j)
levels1 = linspace(-0.000002,0.000002,20);
contourf(X_above,Y_above,log10(eps),([-7 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1,linspace(-0.1,0,5),linspace(0,2,5)]))
titlename = ['|\epsilon| ({\delta}y = ', num2str(y_move(i)*1e6),' {\mu}m, {\delta}x = ', num2str(x_move(i)*1e6),' {\mu}m)'];
title(titlename,'FontSize',24)
xlabel('X (m)','FontSize',18)
ylabel('Y (m)','FontSize',18)
c         = colorbar('FontSize',18);
ylabel(c,'log_{10}(|\epsilon| [{\mu}m])','FontSize',24);
hold on
% % Plot the magnetized material shapes
for k = 1:size(r0,1)
    Plot_magnets(theta0(k),w_vec0(k),h_vec0(k),r0(k,:),1,1,k,N_coil_inf);
end
xlim([-0.1 0.1])
axis equal

if i == 1
figname = 'figure/y + 1um.fig';
elseif i == 2
figname = 'figure/y - 1um.fig';
elseif i == 3
figname = 'figure/x + 1um.fig';
elseif i == 4
figname = 'figure/x - 1um.fig';
end
saveas(gcf,figname)
end

