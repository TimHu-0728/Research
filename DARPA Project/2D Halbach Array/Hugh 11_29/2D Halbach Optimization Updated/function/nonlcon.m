function [c, ceq] = nonlcon(y,target_smooth_length,Mag_width_lower_lim,Mag_width_upper_lim,parameters,N_count,step,M_magnetization,desired_g)
    X          = y(1);              %perturbed length
    Z          = y(2);              %uniform z distance
    lambda     = y(3);              %length of each repeating pattern
    d          = y(4);              %uniform vertical spacing
    z1         = y(5);              %magnet 1 vertical spacing
    z2         = y(6);              %magnet 2 vertical spacing                  
    z3         = y(7);              %magnet 3 vertical spacing
    z4         = y(8);              %magnet 4 vertical spacing
    
    w                   = abs((lambda - 4*d)/4);                                              %mag_width
    h                   = abs(w);                                                             %mag_height
    mag_count           = ceil((target_smooth_length + 2*X + d)/(d+ w));                 %rounding up the number of magnets involved
    mag_overall_length  = mag_count*w + (mag_count-1)*d;
    %% Constructing mags for critical_M
    mags   = zeros(mag_count,6);
    for j = 1:mag_count
        mags(j,1)           = w;
        mags(j,2)           = h;
        mags(j,3)           = -mag_overall_length + (mag_count/2 + 0.5 +(j-1))*w + (mag_count/2 - 0.5+(j-1))*d;
        mags(j:4:end,4)     = 0 - (Z + z1 + h/2);
        mags(j+1:4:end,4)   = 0 - (Z + z2 + h/2);
        mags(j+2:4:end,4)   = 0 - (Z + z3 + h/2);
        mags(j+3:4:end,4)   = 0 - (Z + z4 + h/2);
        mags(j,5)           = pi/2 * (j-1);             % when theta = 0, pointing upwards
        mags(j,6)           = M_magnetization;
    end  
%% preallocating nonlinear inequality constraint and equality constraint
    c          = zeros(1,4+2*N_count);
    ceq        = [];
%% non-linear constraint for the width of magnets  

    c(1)       = w - Mag_width_upper_lim;
    c(2)       = Mag_width_lower_lim - w;
%     c(2)       = 0.003175 - w;

%     c(1)       = w - 1/8* 0.0254*1.01;
%     c(2)       = 1/8*0.0254*0.99 - w;
%     ceq(1)     = w - 1/8*0.0254;
%% non-linear constraint for the number of magnets
%     mag_overall_length = mag_count*w + (mag_count-1)*d;
%     c(3)       = mag_count - floor(0.2/0.003); %upper
%     c(4)       = 0 - mag_count;    %lower limit
    
    c(3)       = ((target_smooth_length + 2*X + d)/(d+ w)) - 200;%floor(0.2/Mag_width_lower_lim); %upper limit of number of mags.
    c(4)       = 0 - ((target_smooth_length + 2*X + d)/(d+ w));                              %lower limit of number of mags.
    
%     c(1)    = mag_overall_length - (target_smooth_length + 2*X);%limits the mag_overall_length to only a particular range
%     c(1)    = mag_overall_length - (target_smooth_length + 2*X + 0.1.);%limits the mag_overall_length to only a particular range
%% non-linear constraint for the 10g requirement
%     [xx, yy, n_x, n_y]                = surface_params(mag_overall_length, N_count, step);
% 
%     [~,~,B_step_up]                   = multiB(xx+ step*n_x,yy+step*n_y,mags);
%     [PI_step_up,~,~,~]     = HB_massForcePotential(xx+ step*n_x , yy+step*n_y,parameters,B_step_up);
% 
%     [~,~,B_step_down]                 = multiB(xx - step*n_x,yy-step*n_y,mags);
%     [PI_step_down,~,~,~] = HB_massForcePotential(xx - step*n_x , yy-step*n_y,parameters,B_step_down);
% 
%     bottom          = 1/(2*step) * ( PI_step_up - PI_step_down);
%     min_bottom_PI   = min(abs(bottom));
%     c(5)            = 10*9.81 - min_bottom_PI; %make it 11g later

%   using an alternative method
    mu0                         = parameters(2);
    [xx, yy, n_x, n_y]          = surface_params(mag_overall_length, N_count, step);
    [~,M,Hx,Hy,~,H_norm]        = critical_M(xx,yy,parameters,n_x,n_y,mags);
    h_dir                       = [Hx;Hy]./H_norm;
    Mx                          = h_dir(1,:).*M;
    My                          = h_dir(2,:).*M;

    [B_r_step_up_wrt_x,B_z_step_up_wrt_x,~]   = multiB(xx + step,yy, mags);
    B_step_up_wrt_x                           = [B_r_step_up_wrt_x;B_z_step_up_wrt_x];

    [B_r_step_down_wrt_x,B_z_step_down_wrt_x,~] = multiB(xx - step,yy, mags);
    B_step_down_wrt_x                           = [B_r_step_down_wrt_x;B_z_step_down_wrt_x];

    [B_r_step_up_wrt_y,B_z_step_up_wrt_y,~]     = multiB(xx, yy+step, mags);
    B_step_up_wrt_y                             = [B_r_step_up_wrt_y;B_z_step_up_wrt_y];

    [B_r_step_down_wrt_y,B_z_step_down_wrt_y,~] = multiB(xx, yy-step, mags);
    B_step_down_wrt_y                           = [B_r_step_down_wrt_y;B_z_step_down_wrt_y];

    H_step_up_wrt_r      = B_step_up_wrt_x/mu0;
    H_step_down_wrt_r    = B_step_down_wrt_x/mu0;
    H_step_up_wrt_z      = B_step_up_wrt_y/mu0;
    H_step_down_wrt_z    = B_step_down_wrt_y/mu0;

    dHdx = 1/(2*step) * (H_step_up_wrt_r - H_step_down_wrt_r);       %equal to zero
    dHdy = 1/(2*step) * (H_step_up_wrt_z - H_step_down_wrt_z);

    f_m             = mu0.* (Mx .* dHdx + My .*dHdy);
%     f_m_mag         = sqrt(f_m(1,:).^2 + f_m(2,:).^2);
    f_m_mag         = abs(f_m(2,:));
    f_m_mag         = f_m_mag/parameters(3);

    c(5:4+N_count)            = desired_g*9.81 - f_m_mag;
%% Rosenweig instability
    [xx, yy, n_x, n_y]  = surface_params(mag_overall_length, N_count, step);
    [Mc,M]              = critical_M(xx,yy,parameters,n_x,n_y,mags);            %calculate M_critical and actual M
    c(5+N_count:end)            = M - Mc + 100;                                         %100 serves as buffer value, essentially Mc -100 >= M


end