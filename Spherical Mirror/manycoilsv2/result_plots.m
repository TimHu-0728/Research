function [actual_contour_x_coordinate, actual_contour_y_coordinate,PI_diff,max_error_in_nm] = result_plots(x, r_flat, I, rz, N_count,parameters)
    % This function is used for all of the plotting that came after the optimization, and used to be on main
    %Calculating some basics, for plotting later
%     N_count         = 91;                                                                                 %numbers of sample points per axis

    Jwave       = customcolormap([0, 0.5, 1], [0, 0, 130; 203,203,203;  134,0,0]/255);
    Jlightwave  = customcolormap([0, 0.5, 1], [0, 0, 1; 1,1,1;  1,0,0]);
    Jcividis    = customcolormap([0, 0.5, 1], [0,33,80; 116,116,116; 252, 230, 71]/255);
    Jaurora     = customcolormap([0, 0.617, 0.77 1], [255,255,255; 29,134,133; 115,0,231; 0,0,135]/255); 
    Jmagma      = customcolormap([0, 0.25, 0.5, 1], [4, 0, 15; 118, 0, 130; 252, 53, 90; 245, 255, 186]/255);
    Jwavei      = flipud(Jwave);
    Jlightwavei = flipud(Jlightwave);
    Jcividisi   = flipud(Jcividis);
    Jaurorai    = flipud(Jaurora);
    Jmagmai     = flipud(Jmagma);



    w               = x(1);
    step            = 0.0001;                                                                             %step size
    
    X_ideal         = zeros(1, N_count);                                                                  %x coordinates of PI_ideal, all zeros
    Y_ideal         = zeros(1, N_count);                                                                  %y coordinates of PI_ideal, all zeros
    [xx, yy, nx, ny] = surface_params(r_flat, N_count, step);                                             %finding points to measure potential, and normal vector away from desired surface at those points
    
    top_top         = massForcePotential(w,xx,yy,I,rz,parameters) ...
        - massForcePotential(w,X_ideal,Y_ideal,I,rz,parameters);
    bottom_bottom   = 1/(2*step) * (massForcePotential(w,xx+step*nx,yy+step*ny,I,rz,parameters) ...
        - massForcePotential(w,xx-step*nx,yy-step*ny,I,rz,parameters));
    error           = top_top./bottom_bottom;
    %delta h
    [max_error,i]   = max(abs(error));                                                                    %max delta h and its location of occurrence
    max_error_in_nm = max_error*10^9                                                                      %deviation of what we have over what we want, difference in frequency
    
    
    %Plotting pi
    PI_ideal        = massForcePotential(w,X_ideal,Y_ideal,I,rz,parameters);
    PI_actual       = massForcePotential(w,xx,yy,I,rz,parameters);
    PI_diff         = PI_actual - PI_ideal;
    XXX             = xx;
    
    %% Deviation on MFP plot

    figure(3)
    hold on
    % xlim([0,15])
    % ylim([-2,9])
    axis equal
    
    %plot(XXX, PI_ideal)
    %plot(XXX, PI_actual,'.')
    plot(XXX(1:end-1), PI_diff(1:end-1),'-')                                                              %TOOK OUT the last element as the value went very high
%     title('Deviation on Mass Force Potential (PI)')
    legend('PI diff')%'PI ideal','PI actual','PI diff')
    xlabel('Radial Displacement (r)')
    ylabel('Value of PI')%('Value of PI')
%     hold off
    
    %LGST Standard
    set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
    set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
        'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
    set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
%     set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
    set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
    set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);
    
    %% Physical Location of Mirror surfaces
    %Plotting the delta H
    ideal_contour_x_coordinate  = xx;
    ideal_contour_y_coordinate  = yy;
    actual_contour_x_coordinate = 0 - error.*nx + ideal_contour_x_coordinate;
    actual_contour_y_coordinate = 0 - error.*ny + ideal_contour_y_coordinate;
    
    figure(4)
    hold on
    % xlim([0,15])
    % ylim([-2,9])
    axis equal
    
    plot(xx, yy,'-.')
    plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'o')
    legend('ideal surface','actual surface')
    title('Physical Location of Mirror surfaces')
    xlabel('Radial Displacement (r)')
    ylabel('Vertical Displacement (z)')

    %LGST Standard
    set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
    set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
        'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
    set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
%     set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
    set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
    set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);

    %%  Plot Potential Contour
    N               = 1000;
    Nq              = 15;
    rlim            = max(rz(:,1))+1;
    rlimq           = max(rz(:,1))+1;
    zlim            = min(rz(:,2))-2;
    shift           = 1;
    [R,Z]           = meshgrid(linspace(0,rlim,N),linspace(-zlim+shift,zlim+shift,N));                % Meshpoints for contour plot
    [Rq,Zq]         = meshgrid(linspace(0,rlimq,Nq),linspace(-zlim+shift,zlim+shift,Nq));              % Meshpoints for quiver plot
    
    [PI,PI_m,PI_w,PI_g,Br,Bz]             = multi_Pi(w,R,Z,I,rz,parameters);
    [PI_q,PI_m_q,PI_w_q,PI_g_q,Br_q,Bz_q] = multi_Pi(w,Rq,Zq,I,rz,parameters);
    [DR_m,DZ_m]                           = gradient(PI_m_q);
    [DR_w,DZ_w]                           = gradient(PI_w_q);
    [DR_g,DZ_g]                           = gradient(PI_g_q);
    [DR,DZ]                               = gradient(PI_q);
    
%     figure(5)
%     title('Optimal Solution')
%     subplot(2,2,1)
%     contourf(R,Z,PI_m,linspace(-20,10,15))          % Adjust contour levels accordingly.
%     h = colorbar;
%     h.Title.String = 'Pi';
%     hold on
%     quiver(Rq,Zq,-DR_m,DZ_m,'Color',[0.4940 0.1840 0.5560])
%     
%     plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
%     plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
%     title('Magnetic Force Potential $\Pi_m$ ','FontSize',25,'Interpreter','latex')
%     ylabel('z (m)','FontSize',20)
%     xlabel('r (m)','FontSize',20)
%     % quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
%     legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',FontSize=10)
%     grid on
% 
%     subplot(2,2,2)
%     contourf(R,Z,PI_g,linspace(-150,50,15))          % Adjust contour levels accordingly.
%     colorbar;
%     hold on
%     quiver(Rq,Zq,DR_g,DZ_g,'Color',[0.4940 0.1840 0.5560])
%     plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
%     plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
%     title('Gravity Force Potential $\Pi_g$ ','FontSize',25,'Interpreter','latex')
%     ylabel('z (m)','FontSize',20)
%     xlabel('r (m)','FontSize',20)
%     % quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
%     legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',FontSize=10)
%     grid on
%     
%     
%     subplot(2,2,3)
%     contourf(R,Z,PI_w,linspace(-150,50,15))          % Adjust contour levels accordingly.
%     colorbar;
%     hold on
%     quiver(Rq,Zq,-DR_w,-DZ_w,'Color',[0.4940 0.1840 0.5560])
%     plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
%     plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
%     title('Centrifugal Force Potential $\Pi_{\omega}$ ','FontSize',25,'Interpreter','latex')
%     ylabel('z (m)','FontSize',20)
%     xlabel('r (m)','FontSize',20)
%     % quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
%     legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',FontSize=10)
%     grid on
%     
%     subplot(2,2,4)
%     contourf(R,Z,PI_g+PI_w,linspace(-150,50,15))          % Adjust contour levels accordingly.
%     colorbar;
%     hold on
%     quiver(Rq,Zq,DR_g-DR_w,DZ_g-DZ_w,'Color',[0.4940 0.1840 0.5560])
%     plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
%     plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
%     title('Centrifugal and Gravity Force Potential $\Pi_g + \Pi_{\omega}$ ','FontSize',25,'Interpreter','latex')
%     ylabel('z (m)','FontSize',20)
%     xlabel('r (m)','FontSize',20)
%     % quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
%     legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',FontSize=10)
%     grid on
    

    %% Magnetic Force Potential
    figure(5)
    title('Optimal Solution')
    contourf(R,Z,PI_m,linspace(-20,10,15))          % Adjust contour levels accordingly.
    h = colorbar;
    h.Title.String = 'Pi';
    hold on
    quiver(Rq,Zq,-DR_m,DZ_m,1,'Color',[0.4940 0.1840 0.5560])%'AutoScale','on','AutoScaleFactor',1.5)
    
    plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
    plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
    title('Magnetic Force Potential $\Pi_m$ ','FontSize',35,'Interpreter','latex')
    ylabel('z (m)','FontSize',35)
    xlabel('r (m)','FontSize',35)
    % quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
    %legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',Location='eastoutside',FontSize=15)
    grid on

    %LGST Standard
    set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
    set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 35,...
        'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
    set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
    %     set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
    set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
    set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);

    colormap(Jcividisi)
    % Save .pdf
    set(gcf,'PaperPositionMode', 'auto','Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'Magnetic Force Potential','-dpdf','-r1000');
    %% Gravitational Force Potential
    figure(6)
    contourf(R,Z,PI_g,linspace(-150,50,15))          % Adjust contour levels accordingly.
    colorbar;
    hold on
    quiver(Rq,Zq,DR_g,DZ_g,'Color',[0.4940 0.1840 0.5560])
    plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
    plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
    title('Gravitational Force Potential $\Pi_g$ ','FontSize',35,'Interpreter','latex')
    ylabel('z (m)','FontSize',35)
    xlabel('r (m)','FontSize',35)
    % quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
    %legend("Force Potential Field Contour",'Force Field','ideal surface','actual surface',Location='eastoutside',FontSize=35)
    grid on
    
    %LGST Standard
    set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
    set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 35,...
        'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
    set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
    %     set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
    set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
    set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);
    colormap(Jcividisi)
    %% Centrifugal Force Potential
    figure(7)
    contourf(R,Z,PI_w,linspace(-150,50,15))          % Adjust contour levels accordingly.
    colorbar;
    hold on
    quiver(Rq,Zq,-DR_w,-DZ_w,'Color',[0.4940 0.1840 0.5560])
    plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
    plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
    title('Centrifugal Force Potential $\Pi_{\omega}$ ','FontSize',35,'Interpreter','latex')
    ylabel('z (m)','FontSize',35)
    xlabel('r (m)','FontSize',35)
    % quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
    %legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',Location='eastoutside',FontSize=15)
    grid on
    %LGST Standard
    set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
    set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 35,...
        'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
    set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
    %     set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
    set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
    set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);
    colormap(Jcividisi)
    %% Centrifugal + Gravitational Force Potential Plot
    figure(8)
    contourf(R,Z,PI_g+PI_w,linspace(-150,50,15))          % Adjust contour levels accordingly.
    colorbar;
    hold on
    quiver(Rq,Zq,DR_g-DR_w,DZ_g-DZ_w,'Color',[0.4940 0.1840 0.5560])
    plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
    plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
    title('Centrifugal and Gravity Force Potential $\Pi_g + \Pi_{\omega}$ ','FontSize',35,'Interpreter','latex')
    ylabel('z (m)','FontSize',35)
    xlabel('r (m)','FontSize',35)
    % quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
    %legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',Location='eastoutside',FontSize=15)
    grid on
    %LGST Standard
    set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
    set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 35,...
        'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
    set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
    %     set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
    set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
    set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);
    colormap(Jcividisi)
    %% Total Force Potential Plot
    figure(9)
    contourf(R,Z,PI,linspace(-100,30,15))          % Adjust contour levels accordingly.
    colorbar;
    hold on
    quiver(Rq,Zq,-DR,DZ,'Color',[0.4940 0.1840 0.5560])
    plot(xx, yy,'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
    plot(actual_contour_x_coordinate, actual_contour_y_coordinate,'bo','LineWidth',1)
    title('Total Force Potential $\Pi$ ','FontSize',35,'Interpreter','latex')
    ylabel('z (m)','FontSize',35)
    xlabel('r (m)','FontSize',35)
    % quiver(Rq,Zq,Br_q,Bz_q,0,'Color',[0.4940 0.1840 0.5560])
    %legend('Force Potential Field Contour','Force Field','ideal surface','actual surface',Location='eastoutside',FontSize=15)
    grid on
    
    %LGST Standard
    set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
    set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 35,...
        'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
    set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
    %     set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
    set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
    set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);
    colormap(Jcividisi)
    %% Plot Magnetization (M) and Critical Magetization (Mc) and Saturation Magnetization (Ms) along the suraface
    figure(10)
    [Mc,M]   = Critical_M(xx',yy',nx',ny',I,rz,parameters);
    plot(xx',Mc,'*',xx',M,'+',xx',ones(size(xx'))*parameters(4),'-','LineWidth',2,'MarkerSize',10)
    title('Magenatization along the interface','FontSize',20)
    xlabel('r (m)','FontSize',18)
    ylabel('M (A/m)','FontSize',18)
    legend('Mc','M','Ms','Location','eastoutside',fontsize=18)
    grid on
        
    %LGST Standard
    set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
    set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
        'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
    set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
    %     set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
    set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
    set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);
    
    % hold off
    %% Calculate magnetic Bond number
%     mui = parameters(2);
%     H
% 
%     Bo_m = mui
end

