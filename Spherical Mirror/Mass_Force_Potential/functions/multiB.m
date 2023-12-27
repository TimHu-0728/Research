function [Br, Bz] = multiB(rold, zold, I, co)%Calculates magnetic field for multiple coild
    %Here, each row of the matrix for the coils used will be of the form [xm ym Bt], using the nomenclature from the B function.
    Br = 0; Bz = 0;
    for i = 1:size(co, 1)
        [Brnew, Bznew] = B(I(i), rold, zold, [co(i, 1) co(i, 2)]);
        Br = Br + Brnew;
        Bz = Bz + Bznew;
    end
end
