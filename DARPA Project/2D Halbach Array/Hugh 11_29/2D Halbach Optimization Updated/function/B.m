function [Bfield_x, Bfield_y] = B(Bt, x_old, y_old, sizee, cenm, t)%Calculates magnetic field for a single magnet
    %x, y is the position in the old coordinate system.
    %sizee = [xe ye] are width and height of magnet, respectively. cenm = [xm, ym] are the coordinates of the magnet's center.
    %Next, we convert both to the "med" coordinate system, which is rotated w.r.t. the original coordinate system by t, the rotation

    medx        = x_old * cos(t) + y_old * sin(t);                    %x pos (of the whole sim. range) (rotated coordinates)
    medy        = y_old * cos(t) - x_old * sin(t);                    %y pos (of the whole sim. range) (rotated coordinates)
    medcenm     = [cenm(1) * cos(t) + cenm(2) * sin(t), cenm(2) * cos(t) - cenm(1) * sin(t)]; %magnet center (rotated coordinates)
    %next, we translate to get the magnet coordinate system
    newx        = medx - medcenm(1);                                             
    newy        = medy - medcenm(2) - sizee(2) / 2;
    lrpp        = 0.5 * log((newx + sizee(1) / 2) .^ 2 + newy .^ 2);             %ln(r) at x + xe / 2, y;
    lrmp        = 0.5 * log((newx - sizee(1) / 2) .^ 2 + newy .^ 2);             %ln(r) at x - xe / 2, y;
    lrpm        = 0.5 * log((newx + sizee(1) / 2) .^ 2 + (newy + sizee(2)) .^ 2);%ln(r) at x + xe / 2, y + ye;
    lrmm        = 0.5 * log((newx - sizee(1) / 2) .^ 2 + (newy + sizee(2)) .^ 2);%ln(r) at x - xe / 2, y + ye;
    dx          = lrpp - lrmp - lrpm + lrmm;                                     %Bx in new coords
    dy          = atan(newy ./ (newx - sizee(1) / 2)) - atan((newy + sizee(2)) ./ (newx - sizee(1) / 2)) - atan(newy ./ (newx + sizee(1) / 2)) + atan((newy + sizee(2)) ./ (newx + sizee(1) / 2));%By in new coords
    Bfield_x     = (dx * cos(t) - dy * sin(t)) * Bt / (2 * pi);                   %B-field in normal coords (rotate it back)
    Bfield_y     = (dy * cos(t) + dx * sin(t)) * Bt / (2 * pi);                   %B-field in normal coords (rotate it back)
    
    %finally, we convert the magnetic field inside the magnet to NaN so it does not show up on the graph.
%     M           = max(0, sign(newx + sizee(1) / 2)) .* max(0, sign(sizee(1) / 2 - newx)) .* max(0, sign(0 - newy)) .* max(0, sign(newy + sizee(2)));
%     M           = 0;                                        %just make them 0
%     Bfieldx = Bfieldx + M;
%     Bfieldy = Bfieldy + M;
end

