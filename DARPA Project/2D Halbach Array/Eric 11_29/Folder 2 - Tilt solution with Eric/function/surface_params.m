function [xx, yy, n_x, n_y] = surface_params(mirror_length, N_count, step)
    %This function generates the points and normal vectors for a surface y = surf_func(x).
    %[IMPORTANT]: surf_func(x) is a parameter of the model and can be changed at will to get a different geometry!!!
%     xx     = linspace(0, mirror_length, N_count);
    xx     = linspace(-mirror_length/2, mirror_length/2, N_count);
    yy     = surf_func(xx);
    
    %generating normal vectors from hat n = hat <del (y - f(x))>
    grad_x = (surf_func(xx + step) - surf_func(xx - step)) / (-2*step);%This is the x-component of the gradient. y-component is 1, so it does not have to be calculated.
    n_x     = grad_x ./ sqrt(grad_x .^ 2 + 1);
    n_y     = 1 ./ sqrt(grad_x .^ 2 + 1);
end