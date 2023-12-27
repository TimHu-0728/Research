function [w,h,mag_count, mag_overall_length, results,mags ] = result_plots(x, parameters, N_count,step,M_magnetization,target_smooth_length, theta2)
    % This function is used for plotting theta
    theta = theta2 * pi / 180;
    results_max = theta * 0;
    results_avg = theta * 0;
    for i = 1:length(theta)
        [w,h,mag_count, mag_overall_length, max_truncated_error_in_nm,average_truncated_error,mags ] = error2(x, parameters, N_count,step,M_magnetization,target_smooth_length, theta(i));
        results_avg(i)=average_truncated_error;
        results_max(i)=max_truncated_error_in_nm;
    end
    results = [results_max, results_avg];

    figure(1);%plotting
    plot(theta2,results_max);
    title('Maximum error vs tilt');
    xlabel('Tilt angle (degrees)');
    ylabel('Maximum error (nm)');

    figure(2);
    plot(theta2,results_avg);
    title('Average error vs tilt');
    xlabel('Tilt angle (degrees)');
    ylabel('Average error (nm)');
end

