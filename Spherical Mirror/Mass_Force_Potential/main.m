%% Initialize parameters and constants
clear
clc
close all

% addpath('functions/')

% Parameters using EMG 700 from Ferrotec corp. Website
% URL: https://ferrofluid.ferrotec.com/products/ferrofluid-emg/water/emg-700-sp/
g     = 9.806;                      % Gravitational Acceleration [m/s^2]
mu0   = 4*pi*1e-7;                  % Vaccum Pemeability [N/A^2]
rho   = 1290;                       % Ferrofluid material density [kg/m^3]
Md    = 28250;                      % Saturation moment of the bulk Ferrofluid [A/m]
d     = 10e-9;                      % Nominal Particle diameter of ferrofluid [m]
phi   = 1;                          % Volume fraction of Ferrofluid
k     = 1.380649e-23;               % Boltzman Constant [m^2*kg/s^2/K]
T     = 293;                        % Temperature [K]
Ms    = 28250;                      % Saturation magnetization of Ferrofluid mateirial [A/m]          
chi0  = 1;                          % Initial Magenetic susceptibility of Ferrofluid (EMG 700 with dilution, SI unit)
gamma = 3*chi0/Ms; 

parameters = [g;mu0;rho;Ms;gamma];



N           = 4;                                                                %number of coils
r           = 30;                                                               %radius of the mirror

%x           = [w , ln(I1) , z1 , r1 , ln(I2) , z2 , r2, ln(I3) , z3 , r3 , ln(I4) , z4 , r4 ]
x0          = [0.5718 , 18 , 0 , 30 , 11  , -1 , 11, 14 , -3 , 20  , 14 , -2.5,25 ];  %initial guess [Hugh]
%enter this in to figure out currents in console after a run: exp([x(2) x(5) x(8) x(11)])
%x           = [1 , 2  ,  3 ,  4 , 5  , 6  , 7 ,  8 ,  9 , 10 , 11  , 12, 13
% x0          = [sqrt(g/r) , 587, -1.5,1.4, 5  , -3 , 13, -200,-1 , 20 , 2000,-2.5,25 ];  %Inspired by Borra

r_upper_lim = 30;
r_lower_lim = 5;    

% z_upper_lim = 5; %%%
% z_lower_lim = 0; %%%
z_upper_lim = 0; 
z_lower_lim = -5;

I_upper_lim = 18;
I_lower_lim = 10;

w_upper_lim = 1*sqrt(g/r);
w_lower_lim = 1*sqrt(g/r);
% w_upper_lim = 5;
% w_lower_lim = 0;

lb          = [w_lower_lim,I_lower_lim,z_lower_lim ,r_lower_lim,I_lower_lim,z_lower_lim ,r_lower_lim,I_lower_lim,z_lower_lim ,r_lower_lim,I_lower_lim,z_lower_lim ,r_lower_lim];
ub          = [w_upper_lim,I_upper_lim,z_upper_lim ,r_upper_lim,I_upper_lim,z_upper_lim ,r_upper_lim,I_upper_lim,z_upper_lim ,r_upper_lim,I_upper_lim,z_upper_lim ,r_upper_lim];

A           = [];                                                                                     %Linear Inequality Constraint, LHS
B           = [];                                                                                     %Linear Inequality Constraint, RHS
Aeq         = [];                                                                                     %Linear equality Constraint, LHS    
Beq         = [];                                                                                     %Linear equality Constraint, RHS


% Options [User Defined]
test     = false;


%% Optimizations and plotting the results
if test
    objective(x0)
    x = [x0(1); exp(x0(2)); x0(3); x0(4); exp(x0(5)); x0(6); x0(7); exp(x0(8)); x0(9); x0(10); exp(x0(11)); x0(12); x0(13)];
    rr = linspace(0, 15, 2e2)
    zz = linspace(-5, 7, 2e2)
    PlotPi(x(1),rr,zz,[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)])
else
    options = optimoptions('fmincon','Algorithm','sqp','FunctionTolerance',1e-12, 'ConstraintTolerance',1e-12,'StepTolerance',1e-12,'OptimalityTolerance',1e-12)
    [x,fval,ef,output,lambda] = fmincon(@objective, x0, A, B, Aeq, Beq, lb, ub, @nlcon,options);          %Starting fmincon 
    [c, ceq]                  = nlcon(x);                                                                 %Calculating our nonlinear inequality and equality constraints
    
    %exponentiating the currents
    x = [x(1); exp(x(2)); x(3); x(4); exp(x(5)); x(6); x(7); exp(x(8)); x(9); x(10); exp(x(11)); x(12); x(13)];

    %Things to check
    disp(x)                                                                                               %Displaying the chosen optimizied values
    sum_RI                    = sum(abs(x(2)*x(4)) + abs(x(5)*x(7)) + abs(x(8)*x(10)) + abs(x(11)*x(13))) %Limit is +-1500000, I*r basically meaning mass of the superconductor
    weighted_PI_difference    = fval                                                                      %user defined objective function value after subbing in x
        
    %Calculating, for plotting 
    N_count         = 91;                                                                                 %numbers of sample points per axis
    t               = linspace(3*pi/2, 3*pi/2 + pi/6, N_count);                                           %angles of bottom right arc of a circle
    xx              = r * cos(t);                                                                         %x coordinates of the simulated range
    yy              = r * sin(t) + r;                                                                     %y coordinates of the simulated range
    X_ideal         = zeros(1, N_count);                                                                  %x coordinates of PI_ideal, all zeros
    Y_ideal         = zeros(1, N_count);                                                                  %y coordinates of PI_ideal, all zeros
    
    top_top         = massForcePotential(x(1),xx,yy,[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters) ...
            - massForcePotential(x(1),X_ideal,Y_ideal,[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters);
    step            = 0.0001;                                                                             %step size
    theta           = atan2d((r-yy),xx);                                                                  %angles used to calculate the normal vector (from +ve x axis CW to the normal vector)
    bottom_bottom   = 1/(2*step) * (massForcePotential(x(1),xx-step*cosd(theta),yy+step*sind(theta),[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters) ...
            - massForcePotential(x(1),xx+step*cosd(theta),yy-step*sind(theta),[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters));
    error           = top_top./bottom_bottom;
    %delta h
    [max_error,i]   = max(abs(error));                                                                         %max delta h and its location of occurrence
    max_error_in_nm = max_error*10^9                                                                      %deviation of what we have over what we want, difference in frequency 
    
    
    %Plotting pi
    PI_ideal        = massForcePotential(x(1),X_ideal,Y_ideal,[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters);
    PI_actual       = massForcePotential(x(1),xx,yy,[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters);
    PI_diff         = PI_actual - PI_ideal;
    XXX             = xx;
    
    
    figure(1)
    hold on
    % xlim([0,15])
    % ylim([-2,9])
    axis equal
    
    %plot(XXX, PI_ideal)
    %plot(XXX, PI_actual,'.')
    plot(XXX, PI_diff,'-')
    title('Deviation on Mass Force Potential (PI)')
    legend('PI diff')%'PI ideal','PI actual','PI diff')
    xlabel('Radial Displacement (r)')
    ylabel('Value of PI')%('Value of PI')
    hold off
    grid on
    
    %Plotting the delta H
    ideal_contour_x_coordinate  = xx;
    ideal_contour_y_coordinate  = yy;
    actual_contour_x_coordinate = error.*cosd(theta) + ideal_contour_x_coordinate;
    actual_contour_y_coordinate = error.*-sind(theta) + ideal_contour_y_coordinate;
    
    figure(2)
    plot(xx, yy,'Color',[0 0.4470 0.7410],'LineWidth',2)
    hold on
    plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'o','LineWidth',1)
    xlim([0,18])
    ylim([-3,5])
    legend('Ideal surface','Actual surface')
    title('Physical Location of Mirror surfaces')
    xlabel('Radial Displacement (r)')
    ylabel('Vertical Displacement (z)')
    grid on
end
% %LGST Standard format
% set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
% set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
%     'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
% set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
% set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
% set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
% set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);

%% Plotting Potentials

r_magnets = [x(4) x(7) x(10) x(13)];
z_magnets = [x(3) x(6) x(9) x(12)];
w         = x(1);

co    = [r_magnets' z_magnets'];
I     = [x(2) x(5) x(8) x(11)]';

N               = 1000;
Nq              = 20;
rlim            = 30;
rlimq           = 20;
zlim            = 6;
shift           = 2;
[R,Z]           = meshgrid(linspace(0,rlim,N),linspace(-zlim+shift,zlim+shift,N));                % Meshpoints for contour plot
[Rq,Zq]         = meshgrid(linspace(0,rlimq,Nq),linspace(-zlim+shift,zlim+shift,Nq));              % Meshpoints for quiver plot

[PI,PI_m,PI_w,PI_g,Br,Bz]             = massForcePotential(w,R,Z,I,co,parameters);
[PI_q,PI_m_q,PI_w_q,PI_g_q,Br_q,Bz_q] = massForcePotential(w,Rq,Zq,I,co,parameters);
[DR_m,DZ_m]                           = gradient(PI_m_q);
[DR_w,DZ_w]                           = gradient(PI_w_q);
[DR_g,DZ_g]                           = gradient(PI_g_q);
[DR,DZ]                               = gradient(PI_q);

% Plotting
figure(3)
title('Optimal Solution')
subplot(2,2,1)
contourf(R,Z,PI_m,linspace(-40,20,15))          % Adjust contour levels accordingly.
h = colorbar;
h.Title.String = 'Pi';
hold on
quiver(Rq,Zq,-DR_m,-DZ_m,'Color',[0.4940 0.1840 0.5560])
plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
title('Magnetic Force Potential $\Pi_m$ ','FontSize',25,'Interpreter','latex')
ylabel('z (m)','FontSize',20)
xlabel('r (m)','FontSize',20)
% quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',FontSize=10)
grid on


subplot(2,2,2)
contourf(R,Z,PI_g,linspace(-150,50,15))          % Adjust contour levels accordingly.
h = colorbar;
h.Title.String = 'Pi';
hold on
quiver(Rq,Zq,-DR_g,-DZ_g,'Color',[0.4940 0.1840 0.5560])
plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
title('Gravity Force Potential $\Pi_g$ ','FontSize',25,'Interpreter','latex')
ylabel('z (m)','FontSize',20)
xlabel('r (m)','FontSize',20)
% quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',FontSize=10)
grid on


subplot(2,2,3)
contourf(R,Z,PI_w,linspace(-150,50,15))          % Adjust contour levels accordingly.
h = colorbar;
h.Title.String = 'Pi';
hold on
quiver(Rq,Zq,-DR_w,-DZ_w,'Color',[0.4940 0.1840 0.5560])
plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
title('Centrifugal Force Potential $\Pi_{\omega}$ ','FontSize',25,'Interpreter','latex')
ylabel('z (m)','FontSize',20)
xlabel('r (m)','FontSize',20)
% quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',FontSize=10)
grid on

subplot(2,2,4)
contourf(R,Z,PI_g+PI_w,linspace(-150,50,15))          % Adjust contour levels accordingly.
h = colorbar;
h.Title.String = 'Pi';
hold on
quiver(Rq,Zq,-DR_g-DR_w,-DZ_g-DZ_w,'Color',[0.4940 0.1840 0.5560])
plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
title('Centrifugal and Gravity Force Potential $\Pi_g + \Pi_{\omega}$ ','FontSize',25,'Interpreter','latex')
ylabel('z (m)','FontSize',20)
xlabel('r (m)','FontSize',20)
% quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',FontSize=10)
grid on


figure(4)
contourf(R,Z,PI,linspace(-150,50,15))          % Adjust contour levels accordingly.
h = colorbar;
h.Title.String = 'Pi';
hold on
quiver(Rq,Zq,-DR,-DZ,'Color',[0.4940 0.1840 0.5560])
plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
title('Total Force Potential $\Pi$ ','FontSize',25,'Interpreter','latex')
ylabel('z (m)','FontSize',20)
xlabel('r (m)','FontSize',20)
% quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',Location='eastoutside',FontSize=15)
grid on



%% Functions

function [c, ceq] = nlcon(x)
    N       = 4;                                                                                       %Number of coils
    % c(1)    = -(abs(x(2)*x(4)) + abs(x(5)*x(7)) + abs(x(8)*x(10)) + abs(x(11)*x(13)));
    % c(2)    = (abs(x(2)*x(4))  + abs(x(5)*x(7)) + abs(x(8)*x(10)) + abs(x(11)*x(13))) - 50*1000*30;
    c       = [];
    ceq     = [];                                                                                      %nonlinear equality constraint
end

function obj = objective(y)
    % Parameters of ferrofluid and constants
    g     = 9.806;                       % Gravitational Acceleration [m/s^2]
    mu0   = 4*pi*1e-7;                   % Vaccum Pemeability [N/A^2]
    rho   = 1290;                        % Ferrofluid material density [kg/m^3]
    Md    = 28250;                      % Saturation moment of the bulk Ferrofluid [A/m]
    d     = 10e-9;                       % Nominal Particle diameter of ferrofluid [m]
    phi   = 1;                           % Volume fraction of Ferrofluid
    k     = 1.380649e-23;                % Boltzman Constant [m^2*kg/s^2/K]
    T     = 293;                         % Temperature [K]
    Ms    = 28250;                      % Saturation magnetization of Ferrofluid mateirial [A/m]          
    gamma = pi/6*mu0*Md*d^3/k/T; 
    chi   = 0.1;                         % Magenetic susceptibility of Ferrofluid
    w     = 0.72;                        % Angular velocity of the container   [rad/s]
    
    parameters = [g;mu0;rho;Ms;gamma];

    x = [y(1); exp(y(2)); y(3); y(4); exp(y(5)); y(6); y(7); exp(y(8)); y(9); y(10); exp(y(11)); y(12); y(13)];
    N_count      = 91;
    r            = 30;                                                                                 %radius
    t            = linspace(3*pi/2, 3*pi/2 + pi/6, N_count);                                           %bottom right arc of circle
    xx           = r * cos(t);
    yy           = r * sin(t) + r;
    X_ideal      = zeros(1, N_count);
    Y_ideal      = zeros(1, N_count);
    top          = massForcePotential(x(1),xx,yy,[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters) ...
        - massForcePotential(x(1),X_ideal,Y_ideal,[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters);
    theta        = atan2d((r-yy),xx);                                                                  %should be from 90:60 degs
    step         = 0.0001;
    bottom       = 1/(2*step) * (massForcePotential(x(1),xx-step*cosd(theta),yy+step*sind(theta),[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters) ...
        - massForcePotential(x(1),xx+step*cosd(theta),yy-step*sind(theta),[x(2);x(5);x(8);x(11)],[x(4) x(3); x(7) x(6); x(10) x(9); x(13) x(12)],parameters));
    error        = top./bottom;
    L2norm       = sqrt(sum(sum(error.^2)));                                                           %L2norm of the error (magnitude)
    Linf         = max(max(abs(error)));                                                               %to account for sudden spikes
    gamma        = 1+4e-10;                                                                                  %tuning parameter for Linf term
    obj          = L2norm/sqrt(max(size(error))) + gamma *Linf;
%     obj          = L2norm + gamma *Linf;
end
