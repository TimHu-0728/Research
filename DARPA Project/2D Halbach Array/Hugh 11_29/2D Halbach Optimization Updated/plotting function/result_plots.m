function [w,h,mag_count, mag_overall_length, max_without_truncated_error_in_nm,average_without_truncated_error, max_truncated_error_in_nm,average_truncated_error,mags ] = result_plots(x, parameters, N_count,step,M_magnetization,plot_x_range, plot_y_range, plot_mag_resolution, plot_quiver_resolution,target_smooth_length)%x, r_flat, I, rz, N_count,step)
    % This function is used for all of the plotting that came after the
    % optimization, as well as calculating all values of interest to
    % evaluate the performance of the optimized design

    X          = x(1);              %perturbed length
    Z          = x(2);              %uniform z distance
    lambda     = x(3);              %length of each repeating pattern
    d          = x(4);              %uniform vertical spacing
    z1         = x(5);              %magnet 1 vertical spacing
    z2         = x(6);              %magnet 2 vertical spacing                  
    z3         = x(7);              %magnet 3 vertical spacing
    z4         = x(8);              %magnet 4 vertical spacing

    w                  = (lambda - 4*d)/4;                                    %mag_width
    h                  = w;                                                   %mag_height, assuming its a square magnet for now
    mag_count          = ceil((target_smooth_length + 2*X + d)/(d+ w));       %rounding up the number of magnets involved
    mag_overall_length = mag_count*w + (mag_count-1)*d;

%% Construction of mags
    % mags = [uniform width, uniform height, x_center, y_center, theta, magnetization]
    mags  = zeros(mag_count,6);
    for j = 1:mag_count
        mags(j,1)           = w;
        mags(j,2)           = h;
        mags(j,3)           = -mag_overall_length + (mag_count/2 + 0.5 +(j-1))*w + (mag_count/2 - 0.5+(j-1))*d;
        mags(j:4:end,4)     = 0 - (Z + z1 + h/2);
        mags(j+1:4:end,4)   = 0 - (Z + z2 + h/2);
        mags(j+2:4:end,4)   = 0 - (Z + z3 + h/2);
        mags(j+3:4:end,4)   = 0 - (Z + z4 + h/2);
        mags(j,5)           = pi/2 * (j-1);                                   % when theta = 0, pointing upwards
        mags(j,6)           = M_magnetization;
    end  
    
    X_ideal             = zeros(1, N_count);
    Y_ideal             = zeros(1, N_count);
    [xx, yy, n_x, n_y]  = surface_params(mag_overall_length, N_count, step);
%     [xx, yy, n_x, n_y]  = surface_params(mirror_length, N_count, step);

    %% Plotting pi
    [~,~,B_ideal]       = multiB(X_ideal,Y_ideal,mags);
    PI_ideal            = HB_massForcePotential(X_ideal, Y_ideal,parameters,B_ideal);

    [~,~,B_actual]      = multiB(xx,yy,mags);
    PI_actual           = HB_massForcePotential(xx,yy,parameters,B_actual);

    PI_diff             = PI_actual - PI_ideal;
    top                 = PI_diff;

    [~,~,B_step_up]     = multiB(xx + step*n_x,yy+step*n_y,mags);
    PI_step_up          = HB_massForcePotential(xx+ step*n_x , yy+step*n_y,parameters,B_step_up);

    [~,~,B_step_down]   = multiB(xx - step*n_x,yy-step*n_y,mags);
    PI_step_down        = HB_massForcePotential(xx - step*n_x , yy-step*n_y,parameters,B_step_down);

    bottom              = 1/(2*step) * (PI_step_up - PI_step_down);
    error               = top./bottom;


    %delta h
    [max_error,i]                     = max(abs(error));                                    %max delta h and its location of occurrence
    max_without_truncated_error_in_nm = max_error*10^9 ;                                    %deviation of what we have over what we want, difference in nm
    average_without_truncated_error   = mean(abs(error)) * 10^9   ;                         %average error, in nm

    %% Truncated error
%     [xx, yy, n_x, n_y]  = surface_params(mag_overall_length-2*X, N_count, step);            %The only difference, as now the overall length has been reduced by 2X
% 
%     [~,~,B_ideal]       = multiB(X_ideal,Y_ideal,mags);
%     PI_ideal            = HB_massForcePotential(X_ideal, Y_ideal,parameters,B_ideal);
% 
%     [~,~,B_actual]      = multiB(xx,yy,mags);
%     PI_actual           = HB_massForcePotential(xx,yy,parameters,B_actual);
% 
%     PI_diff             = PI_actual - PI_ideal;
%     top                 = PI_diff;
% 
%     [~,~,B_step_up]     = multiB(xx + step*n_x,yy+step*n_y,mags);
%     PI_step_up          = HB_massForcePotential(xx+ step*n_x , yy+step*n_y,parameters,B_step_up);
% 
%     [~,~,B_step_down]   = multiB(xx - step*n_x,yy-step*n_y,mags);
%     PI_step_down        = HB_massForcePotential(xx - step*n_x , yy-step*n_y,parameters,B_step_down);
% 
%     bottom              = 1/(2*step) * (PI_step_up - PI_step_down );
%     error               = top./bottom;
    
%     starting_index  = round(N_count/mag_overall_length * X);
    starting_index  = ceil(N_count/mag_overall_length * X);
    ending_index    = round(N_count/mag_overall_length * (mag_overall_length - X));
%     ending_index    = floor(N_count/mag_overall_length * (mag_overall_length - X));
    truncated_error = error(starting_index:ending_index);
    [max_truncated_error,i]   = max(abs(truncated_error));                                                                         %max delta h and its location of occurrence
    max_truncated_error_in_nm = max_truncated_error*10^9 ;                                                                          %deviation of what we have over what we want, difference in nm
    average_truncated_error   = mean(abs(truncated_error)) * 10^9  ;                                                                %average error, in nm
    

    %% Plotting
    
    %For deltaPI
    figure(1)
    hold on
    axis equal
    plot(xx, PI_diff,'-')
    title('Deviation on Mass Force Potential (PI)')
    legend('PI diff')                                       %'PI ideal','PI actual','PI diff')
    xlabel('X Displacement (m)')
    ylabel('Value of PI')                                   %('Value of PI')
    hold off
    
    %For actual surface
    %Plotting the delta H
    ideal_contour_x_coordinate  = xx;
    ideal_contour_y_coordinate  = yy;
    actual_contour_x_coordinate = 0 - error.*n_x + ideal_contour_x_coordinate;
    actual_contour_y_coordinate = 0 - error.*n_y + ideal_contour_y_coordinate;
    
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
    
    %For halbach array
    figure(3)
    plotB(-plot_x_range, -plot_y_range, plot_x_range, plot_y_range, plot_mag_resolution, plot_quiver_resolution, mags);
end

