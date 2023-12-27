clc
clear
close all
addpath('./used_function','./figure')

%%
% Initialize parameters  [User Input]
tic
lambda      = 0.25*0.0254; [0.5 3/2 2 4]*0.0254;
w           = lambda/4;                                     % [m]     Width of magnet
h           = lambda/2;                                     % [m]     Height of magnet 
desired_len = [0.5 0.25 0.15];                              % [m]     Desired length of the halbach array
gridsz      = 3000;                                         % Grid size for contour plot, increase this to get more accuracy but need more computation cost
h_ff2hb     = 0.004;                                        % [m]     Ferrofluid height h
g           = 9.81;                                         % [m/s^2] Gravity 
mu0         = 4*pi*1e-7; 
N_coil      = zeros(3,3);
chi0        = 0.181;
[rho,Ms]    = get_Rho_Ms(chi0);                                         %obtain rho and Ms specs from Chi0 interpolation
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

% User this for simulation field above the magnets 
% y_ind_above = find(rangey <= -h_ff2hb,1,'last');
% Use this if you want the simulation field to include magnets
y_ind_above = 1;
% Set this to 1 if want to plot PI-PI_ref, 0 if want to plot only PI. PI is
% the total potential [m^2/s^2]
minus_PI_ref = 1;

figure(1)
hold on
for j = 1:length(lambda)
    
    if minus_PI_ref
        % Calculate infinite halbach array PI
        N_coil_inf  = 1000;
        Mm          = (1.48/mu0)*ones(N_coil_inf,1);                    % [A/m]   Magnetization
        Ke          = Mm;                                                % [A/m]   Sheet current         

        %magnet placement (x,y coordinates)
        r0                              = zeros(N_coil_inf,2);
        r0(1:N_coil_inf/2,1)            = w(j)/2;
        r0(N_coil_inf/2+1:N_coil_inf,1) = -w(j)/2;
        r0(1:N_coil_inf/2,1)            = r0(1:N_coil_inf/2,1) - w(j)*(N_coil_inf/2:-1:1)';
        r0(N_coil_inf/2+1:N_coil_inf,1) = r0(N_coil_inf/2+1:N_coil_inf,1) + w(j)*(1:1:N_coil_inf/2)';
        r0(:,2)                         = -h(j)/2-h_ff2hb;
        %magnet magnetization direction
        theta           = zeros(N_coil_inf,1);
        theta(1:4:end)  = deg2rad(0);
        theta(2:4:end)  = deg2rad(90);
        theta(3:4:end)  = deg2rad(180);
        theta(4:4:end)  = deg2rad(270);

        %width input for calculating B field
        w_vec           = ones(N_coil_inf,1)*w(j);
        w_vec(2:2:end)  = h(j);

        %height input for calculating B field
        h_vec           = ones(N_coil_inf,1)*h(j);
        h_vec(2:2:end)  = w(j);
        % Calculates B fields [T]
        [Bx,By]         = calculatingB(mu0,Ke,h_vec,w_vec,X,Y,theta,r0);         % B field for contour plot
        B               = sqrt(Bx.^2+By.^2);  
        Hx              = Bx / mu0;                                                 %obtaining X magnetic field from B field
        Hy              = By / mu0;                                                 %obtaining Y magnetic field from B field
        H_norm          = sqrt(Hx.^2+Hy.^2);
        % Plot 2D HB if need if need to investigate

        % PI_ref
        X_above     = X(y_ind_above:end,:);
        Y_above     = Y(y_ind_above:end,:);
        H_norm_above= H_norm(y_ind_above:end,:); 
        rangey_above= rangey(y_ind_above:end);
        theta_0     = deg2rad(0);
        PI_m        = -mu0/rho*Ms*lnsinh(gamma * H_norm_above)/ gamma;
        PI_g        = g*(Y_above.* cos(theta_0) - X_above.* sin(theta_0)); 
        PI_ref      = PI_g+PI_m;
        % Gradient of PI_ref
        [nabla_PI_ref_x,nabla_PI_ref_y] = gradient(PI_ref, hx_eric, hy_eric);    
        levels      = linspace(-1,1,40);            % Set this accordingly   
    else
        PI_ref      = zeros(size(X));
        levels      = linspace(-5,1,30);
    end
    % Plot PI_ref 3 plots
    %magnet placement (x,y coordinates)

        N_coil_inf   = ceil(0.6/w(j));
       if mod(N_coil_inf,2) == 1                                       % Make even number of magnets
        N_coil_inf = N_coil_inf+1;
       end
        r0                              = zeros(N_coil_inf,2);
        r0(1:N_coil_inf/2,1)            = w(j)/2;
        r0(N_coil_inf/2+1:N_coil_inf,1) = -w(j)/2;
        r0(1:N_coil_inf/2,1)            = r0(1:N_coil_inf/2,1) - w(j)*(N_coil_inf/2:-1:1)';
        r0(N_coil_inf/2+1:N_coil_inf,1) = r0(N_coil_inf/2+1:N_coil_inf,1) + w(j)*(1:1:N_coil_inf/2)';
        r0(:,2)                         = -h(j)/2-h_ff2hb;
        %magnet magnetization direction
        theta           = zeros(N_coil_inf,1);
        theta(1:4:end)  = deg2rad(0);
        theta(2:4:end)  = deg2rad(90);
        theta(3:4:end)  = deg2rad(180);
        theta(4:4:end)  = deg2rad(270);

        %width input for calculating B field
        w_vec           = ones(N_coil_inf,1)*w(j);
        w_vec(2:2:end)  = h(j);

        %height input for calculating B field
        h_vec           = ones(N_coil_inf,1)*h(j);
        h_vec(2:2:end)  = w(j);
     
        % Plot PI_inf Contour
        figure(1)
        %subplot(2,2,j)
        levels1 = linspace(-8,0,20);
        contourf(X_above,Y_above,PI_ref,levels1)
        titlename = ['\lambda = ', num2str(lambda(j)/0.0254),' in.'];
        title(titlename,'FontSize',18)
        xlabel('X (m)','FontSize',18)
        ylabel('Y (m)','FontSize',18)
        c              = colorbar('FontSize',18);
        ylabel(c,'PI [m^2/s^2]','FontSize',18);
        hold on
        % % Plot the magnetized material shapes
        for k = 1:size(r0,1)
            Plot_magnets(theta(k),w_vec(k),h_vec(k),r0(k,:),1,0);
        end
        xlim([-0.3 0.3])
        axis equal

    for i = 1:length(desired_len)
        Mm          = (1.48/mu0)*ones(N_coil(i,j),1);                    % [A/m]   Magnetization
        Ke          = Mm;                                                % [A/m]   Sheet current         
        
        %magnet placement (x,y coordinates)
        r0                                = zeros(N_coil(i,j),2);
        r0(1:N_coil(i,j)/2,1)             = w(j)/2;
        r0(N_coil(i,j)/2+1:N_coil(i,j),1) = -w(j)/2;
        r0(1:N_coil(i,j)/2,1)             = r0(1:N_coil(i,j)/2,1) - w(j)*(N_coil(i,j)/2:-1:1)';
        r0(N_coil(i,j)/2+1:N_coil(i,j),1) = r0(N_coil(i,j)/2+1:N_coil(i,j),1) + w(j)*(1:1:N_coil(i,j)/2)';
        r0(:,2)                           = -h(j)/2-h_ff2hb;
        %magnet magnetization direction
        theta           = zeros(N_coil(i,j),1);
        theta(1:4:end)  = deg2rad(0);
        theta(2:4:end)  = deg2rad(90);
        theta(3:4:end)  = deg2rad(180);
        theta(4:4:end)  = deg2rad(270);
        
        %width input for calculating B field
        w_vec           = ones(N_coil(i,j),1)*w(j);
        w_vec(2:2:end)  = h(j);
        
        %height input for calculating B field
        h_vec           = ones(N_coil(i,j),1)*h(j);
        h_vec(2:2:end)  = w(j);

        % Calculates B fields [T]
        [Bx,By]         = calculatingB(mu0,Ke,h_vec,w_vec,X,Y,theta,r0);         % B field for contour plot
        B               = sqrt(Bx.^2+By.^2);  
        Hx              = Bx / mu0;                                                 %obtaining X magnetic field from B field
        Hy              = By / mu0;                                                 %obtaining Y magnetic field from B field
        H_norm          = sqrt(Hx.^2+Hy.^2);
        % Plot 2D HB if need if need to investigate

        % PI
        X_above     = X(y_ind_above:end,:);
        Y_above     = Y(y_ind_above:end,:);
        H_norm_above= H_norm(y_ind_above:end,:); 
        rangey_above= rangey(y_ind_above:end);
        theta_0     = deg2rad(0);
        PI_m        = -mu0/rho*Ms*lnsinh(gamma * H_norm_above)/ gamma;
        PI_g        = g*(Y_above.* cos(theta_0) - X_above.* sin(theta_0)); 
        PI          = PI_g+PI_m-PI_ref;
        deviation_y = (PI)./nabla_PI_ref_y;
        

        figure
        %subplot(2,2,j)
        levels1 = linspace(-8,0,20);
        contourf(X_above,Y_above,PI,levels1)
        titlename = ['\lambda = ', num2str(lambda(j)/0.0254),' in.'];
        title(titlename,'FontSize',18)
        xlabel('X (m)','FontSize',18)
        ylabel('Y (m)','FontSize',18)
        c              = colorbar('FontSize',18);
        ylabel(c,'PI [m^2/s^2]','FontSize',18);
        hold on
        % % Plot the magnetized material shapes
        for k = 1:size(r0,1)
            Plot_magnets(theta(k),w_vec(k),h_vec(k),r0(k,:),1,0);
        end
        xlim([-0.3 0.3])
        axis equal
       
        % Plot PI, deviation
        figure
        %subplot(3,3,j+3*(i-1)) 
        % log10(1e6 * abs(deviation_y)), [-2:1:3]
        contourf(X_above,Y_above,deviation_y)
        title_name     = ['Edge effects (L = ', num2str(desired_len(i)*100),' cm , ','\lambda = ',num2str(lambda(j)/0.0254),' in.)'];
        title(title_name,'FontSize',18)
        if i == 3
            xlabel('X (m)','FontSize',18)
        end
        if j == 1
            ylabel('Y (m)','FontSize',18)
        end
        c              = colorbar('FontSize',18);
        ylabel(c,'log_{10}(\epsilon [\mum])','FontSize',18);
        clim([-2,3])
        hold on

        % Plot the magnetized material shapes
        for k = 1:size(r0,1)
            Plot_magnets(theta(k),w_vec(k),h_vec(k),r0(k,:),1,0);
        end
        axis equal
        
        % Plot format
        set(findobj('Type', 'Legend'), 'Interpreter', 'tex', 'box','on');
        set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
            'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'tex','TickDir','out');
        set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'tex');
        set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [2,2,25,15]);
        set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'tex');
        drawnow

        % Save
        filename1 = ['./figure/','Edge effects (L = ', num2str(desired_len(i)*100),' cm , ','lambda = ',num2str(lambda(j)/0.0254),' in.)','.png'];
        saveas(gcf, filename1);

    end
end
% figure(1)
% sgtitle('Deviation from ideal equipotential line',fontsize = 24)
% filename1 = 'figure/PI_analysis.png'; 
% saveas(gcf, filename1); 
% figure(2)
% sgtitle('PI (inifinitely long)',fontsize = 24)
% filename2 = 'figure/PI_infinite_long.png'; 
% saveas(gcf, filename2); 
toc

