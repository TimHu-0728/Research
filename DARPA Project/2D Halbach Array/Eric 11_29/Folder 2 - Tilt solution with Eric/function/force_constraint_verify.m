function [f_m_mag,violating_index,B_mag] = force_constraint_verify(mag_tot_length,X,N_count,step,parameters,mags,desired_g)
mu0                         = parameters(2);
[xx, yy, n_x, n_y]          = surface_params(mag_tot_length , N_count, step);
% [xx, yy, n_x, n_y]          = surface_params(mag_tot_length - 2*X, N_count, step);
[~,M,Hx,Hy,~,H_norm]        = critical_M(xx,yy,parameters,n_x,n_y,mags);
h_dir                       = [Hx;Hy]./H_norm;
Mx                          = h_dir(1,:).*M;
My                          = h_dir(2,:).*M;

%Calculating dHdx and dHdy
% [B_r_step_up, B_z_step_up,~]     = multiB(xx + step*n_x,yy+step*n_y, mags)
% [B_r_step_down, B_z_step_down,~] = multiB(xx - step*n_x,yy-step*n_y, mags)
% [~,~,B_mag_step_up_wrt_x]   = multiB(xx + step,yy, mags);
% [~,~,B_mag_step_down_wrt_x] = multiB(xx - step,yy, mags);
% [~,~,B_mag_step_up_wrt_y]   = multiB(xx, yy+step, mags);
% [~,~,B_mag_step_down_wrt_y] = multiB(xx, yy-step, mags);


[B_r_step_up_wrt_x,B_z_step_up_wrt_x,~]   = multiB(xx + step,yy, mags);
B_step_up_wrt_x = [B_r_step_up_wrt_x;B_z_step_up_wrt_x];

[B_r_step_down_wrt_x,B_z_step_down_wrt_x,~] = multiB(xx - step,yy, mags);
B_step_down_wrt_x = [B_r_step_down_wrt_x;B_z_step_down_wrt_x];

[B_r_step_up_wrt_y,B_z_step_up_wrt_y,~]   = multiB(xx, yy+step, mags);
B_step_up_wrt_y = [B_r_step_up_wrt_y;B_z_step_up_wrt_y];

[B_r_step_down_wrt_y,B_z_step_down_wrt_y,~] = multiB(xx, yy-step, mags);
B_step_down_wrt_y = [B_r_step_down_wrt_y;B_z_step_down_wrt_y];


H_step_up_wrt_r      = B_step_up_wrt_x/mu0;
H_step_down_wrt_r    = B_step_down_wrt_x/mu0;
H_step_up_wrt_z      = B_step_up_wrt_y/mu0;
H_step_down_wrt_z    = B_step_down_wrt_y/mu0;
% 
dHdx = 1/(2*step) * (H_step_up_wrt_r - H_step_down_wrt_r);       %equal to zero
dHdy = 1/(2*step) * (H_step_up_wrt_z - H_step_down_wrt_z);

f_m = mu0.* (Mx .* dHdx + My .*dHdy);
% f_m_mag = sqrt(f_m(1,:).^2 + f_m(2,:).^2);
f_m_mag         = abs(f_m(2,:));
violating_index = find(abs(f_m_mag) < 10*9.81);

B = mu0 .*([Mx;My] + [Hx;Hy]);
B_mag = sqrt(B(1,:).^2 + B(2,:).^2);

if min(abs(f_m_mag))<desired_g*9.81
    disp('10g constraint violated')
    for i = 1:length(violating_index)
        if xx(violating_index(i)) <= -mag_tot_length/2 + X || xx(violating_index(i)) >= mag_tot_length/2 - X
            disp(["The constraints violation occurs outside range of interest, index:", num2str(violating_index(i))])
        else
            disp(["The constraint violation occurs WITHIN the range of interest, index:", num2str(violating_index(i))])
        end
    end
else
    disp('constraints satisfied')
end

end

