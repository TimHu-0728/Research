function [w,h,mag_count, mag_overall_length, max_truncated_error_in_nm,average_truncated_error,mags ] = error2(x, parameters, N_count,step,M_magnetization,target_smooth_length, theta)
    % This function is used for findig the error of a tilted plate.
    
    X                       = x(1);              %perturbed length
    Z                       = x(2);              %uniform z distance
    lambda                  = x(3);              %length of each repeating pattern
    d                       = x(4);              %uniform vertical spacing
    z1                      = x(5);              %magnet 1 vertical spacing
    z2                      = x(6);              %magnet 2 vertical spacing
    z3                      = x(7);              %magnet 3 vertical spacing
    z4                      = x(8);              %magnet 4 vertical spacing
    
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

    %% Pi
    [~,~,B_ideal]       = multiB(X_ideal,Y_ideal,mags);
    PI_ideal            = HB_massForcePotential(X_ideal, Y_ideal,parameters,B_ideal, theta);

    [~,~,B_actual]      = multiB(xx,yy,mags);
    PI_actual           = HB_massForcePotential(xx,yy,parameters,B_actual, theta);

    PI_diff             = PI_actual - PI_ideal;
    top                 = PI_diff;

    [~,~,B_step_up]     = multiB(xx + step*n_x,yy+step*n_y,mags);
    PI_step_up          = HB_massForcePotential(xx+ step*n_x , yy+step*n_y,parameters,B_step_up, theta);

    [~,~,B_step_down]   = multiB(xx - step*n_x,yy-step*n_y,mags);
    PI_step_down        = HB_massForcePotential(xx - step*n_x , yy-step*n_y,parameters,B_step_down, theta);

    bottom              = 1/(2*step) * (PI_step_up - PI_step_down);
    error               = top./bottom;


    %delta h
    [max_error,i]                     = max(abs(error));                                    %max delta h and its location of occurrence

    starting_index  = ceil(N_count/mag_overall_length * X);
    ending_index    = round(N_count/mag_overall_length * (mag_overall_length - X));
%     ending_index    = floor(N_count/mag_overall_length * (mag_overall_length - X));
    truncated_error = error(starting_index:ending_index);
    [max_truncated_error,i]   = max(abs(truncated_error));                                                                         %max delta h and its location of occurrence
    max_truncated_error_in_nm = max_truncated_error*10^9 ;                                                                          %deviation of what we have over what we want, difference in nm
    average_truncated_error   = mean(abs(truncated_error)) * 10^9;                                                                 %average error, in nm
end

