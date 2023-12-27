function [B_x, B_y,B_mag] = multiB(x_old, y_old, mags)%Calculates magnetic field for multiple magnets
    %Here, each row of the matrix for the magnets used will be of the form [xe ye xm ym t Bt], using the nomenclature from the B function.
    B_x            = 0;
    B_y            = 0;
    for i = 1:size(mags, 1)
        [Bx_new, By_new] = B(mags(i, 6), x_old, y_old, [mags(i, 1) mags(i, 2)], [mags(i, 3) mags(i, 4)], mags(i, 5));
        B_x             = B_x + Bx_new;
        B_y             = B_y + By_new;
    end
        B_mag           = sqrt( B_x.^2 + B_y.^2);
end

