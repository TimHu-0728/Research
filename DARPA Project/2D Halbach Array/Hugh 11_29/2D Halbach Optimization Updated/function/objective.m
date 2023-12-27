function obj = objective(y, parameters, N_count,step,M_magnetization,target_smooth_length,desired_g)
%y = [X, Z, lambda, d, z1, z2, z3, z4]
    X          = y(1);              %perturbed length
    Z          = y(2);              %uniform z distance
    lambda     = y(3);              %length of each repeating pattern
    d          = y(4);              %uniform vertical spacing
    z1         = y(5);              %magnet 1 vertical spacing
    z2         = y(6);              %magnet 2 vertical spacing                  
    z3         = y(7);              %magnet 3 vertical spacing
    z4         = y(8);              %magnet 4 vertical spacing

    w                  = (lambda - 4*d)/4;                                      %mag_width
    h                  = w;                                                     %mag_height, assuming its a square magnet for now
    mag_count          = ceil((target_smooth_length + 2*X + d)/(d+ w));         %rounding up the number of magnets involved
    mag_overall_length = mag_count*w + (mag_count-1)*d;

    %% Constructing mags
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
%         mags(j,3) = (w/2 + (j-1)*w) + (j-1)*d - (w * (4*Cycle_length) + d * (4*Cycle_length - 1))/2 ;
%         mags(j,4) = 0 - (y(4+j) + h/2);
        mags(j,5)           = pi/2 * (j-1);                             % when theta = 0, pointing upwards
        mags(j,6)           = M_magnetization;
    end 
%% optimizing for the min. error over entire region
    [xx, yy, n_x, n_y]  = surface_params(mag_overall_length, N_count, step);

%% optimizing for the min. error over truncated region
%     [xx, yy, n_x, n_y]  = surface_params(mag_overall_length-2*X, N_count, step);

%% Objective function
    X_ideal             = zeros(1, N_count);
    Y_ideal             = zeros(1, N_count);

    [~,~,B_ideal]   = multiB(X_ideal,Y_ideal,mags);
    PI_ideal        = HB_massForcePotential(X_ideal, Y_ideal,parameters,B_ideal);

    [~,~,B_actual]  = multiB(xx,yy,mags);
    PI_actual       = HB_massForcePotential(xx,yy,parameters,B_actual);

    top             = PI_actual - PI_ideal;

    [~,~,B_step_up]                   = multiB(xx+ step*n_x,yy+step*n_y,mags);
    [PI_step_up,PI_step_up_m,~,~]     = HB_massForcePotential(xx+ step*n_x , yy+step*n_y,parameters,B_step_up);

    [~,~,B_step_down]                 = multiB(xx - step*n_x,yy-step*n_y,mags);
    [PI_step_down,PI_step_down_m,~,~] = HB_massForcePotential(xx - step*n_x , yy-step*n_y,parameters,B_step_down);


    bottom       = 1/(2*step) * ( PI_step_up - PI_step_down);
    error        = top./bottom;

    %% additional penalty function for the 10g constraint
%     bottom_m          = 1/(2*step) * ( PI_step_up_m - PI_step_down_m);
%     min_bottom_PI_m   = min(abs(bottom_m));                                                            %choosing out the smallest PI_m available
    min_bottom_PI_m = min(abs(bottom));
    if min_bottom_PI_m >= desired_g*9.81
        acceleration_penalty = 0;                                                                      %10g constraint satisfied
    else
%         acceleration_penalty = exp(abs(100));
%           acceleration_penalty = exp(abs(200));
%         acceleration_penalty = 3.2690e+06;
%         acceleration_penalty = 5e+24;
        acceleration_penalty = 1e-4*(desired_g*9.81-min_bottom_PI_m);
%         acceleration_penalty = exp(abs(150 - min_bottom_PI_m));                                          %10g constraint not satisfied
    end 

    %% additional penalty function for having not correct dim. of magnet

%     discretized_values  = [1/8,3/16,1/4,3/8,1/2,3/4,1,1.5]*0.0254;                                     %0.003175	0.0047625	0.00635	0.009525	0.0127	0.01905	0.0254	0.0381
    discretized_values  = [1/8]*0.0254;
%     if any(abs(w - discretized_values) < 1e-6)                                                         % within the tolerance of discretized mag dimension
%         sizing_penalty  = 0;
%     else
%         sizing_penalty  = 1e+3;
% %         sizing_penalty  = exp(abs(100));
%     end
    sizing_penalty = max(0, 1e3*(0-w)-h*1e3);
   
    L2norm       = sqrt(sum(sum(error.^2)));                                                           %L2norm of the error (magnitude)
    Linf         = max(max(abs(error)));                                                               %to account for sudden spikes
    gamma        = 1;                                                                                  %tuning parameter for Linf term
    acceleration_penalty
    obj          = L2norm/sqrt(max(size(error))) + gamma *Linf + acceleration_penalty + sizing_penalty


%       obj        = acceleration_penalty + sizing_penalty
%     obj          = sizing_penalty 
%     obj = L2norm/sqrt(max(size(error))) + gamma *Linf 
%       obj = L2norm/sqrt(max(size(error))) + gamma *Linf + sizing_penalty;
end