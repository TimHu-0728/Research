function result_plots(x, r_flat, I, rz, N_count)
    % This function is used for all of the plotting that came after the optimization, and used to be on main
    %Calculating some basics, for plotting later
%     N_count         = 91;                                                                                 %numbers of sample points per axis
    
    % Parameters using EMG 700 from Ferrotec corp. Website
    % URL: https://ferrofluid.ferrotec.com/products/ferrofluid-emg/water/emg-700-sp/
    g     = 9.806;                      % Gravitational Acceleration [m/s^2]
    mu0   = 4*pi*1e-7;                  % Vaccum Pemeability [N/A^2]
    rho   = 1290;                       % Ferrofluid material density [kg/m^3]
    Ms    = 28250;                      % Saturation magnetization of Ferrofluid mateirial [A/m]          
    chi0  = 1;                          % Initial Magenetic susceptibility of Ferrofluid (EMG 700 with dilution, SI unit)
    gamma = 3*chi0/Ms; 
    
    parameters = [g;mu0;rho;Ms;gamma];



    step            = 0.0001;                                                                             %step size
    
    X_ideal         = zeros(1, N_count);                                                                  %x coordinates of PI_ideal, all zeros
    Y_ideal         = zeros(1, N_count);                                                                  %y coordinates of PI_ideal, all zeros
    [xx, yy, nx, ny] = surface_params(r_flat, N_count, step);                                             %finding points to measure potential, and normal vector away from desired surface at those points
    
    top_top         = massForcePotential(x(1),xx,yy,I,rz,parameters) ...
        - massForcePotential(x(1),X_ideal,Y_ideal,I,rz,parameters);
    bottom_bottom   = 1/(2*step) * (massForcePotential(x(1),xx+step*nx,yy+step*ny,I,rz,parameters) ...
        - massForcePotential(x(1),xx-step*nx,yy-step*ny,I,rz,parameters));
    error           = top_top./bottom_bottom;
    %delta h
    [max_error,i]   = max(abs(error));                                                                         %max delta h and its location of occurrence
    max_error_in_nm = max_error*10^9                                                                      %deviation of what we have over what we want, difference in frequency
    
    
    %Plotting pi
    PI_ideal        = massForcePotential(x(1),X_ideal,Y_ideal,I,rz,parameters);
    PI_actual       = massForcePotential(x(1),xx,yy,I,rz,parameters);
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
    
    
    %Plotting the delta H
    ideal_contour_x_coordinate  = xx;
    ideal_contour_y_coordinate  = yy;
    actual_contour_x_coordinate = 0 - error.*nx + ideal_contour_x_coordinate;
    actual_contour_y_coordinate = 0 - error.*ny + ideal_contour_y_coordinate;
    
    figure(2)
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
    hold off
end

