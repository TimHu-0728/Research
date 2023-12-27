function plotB(x1, y1, x2, y2, N, Nq, mags)%plots (for now) - q is for quiver
    x        = linspace(x1, x2, N);                   %x coord matrix
    y        = linspace(y1, y2, N);                   %y coord matrix
    pre_xq   = linspace(x1, x2, Nq + 1);              %x coord for quiver matrix
    pre_yq   = linspace(y1, y2, Nq + 1);              %y coord for quiver matrix
    xq       = (pre_xq(1:Nq) + pre_xq(2:Nq + 1)) / 2; %only use the center x coord for quiver
    yq       = (pre_yq(1:Nq) + pre_yq(2:Nq + 1)) / 2; %only use the center y coord for quiver
    [X, Y]   = meshgrid(x, y);                        
    [U, V]   = multiB(X, Y, mags);                    %generating magnitude output (contour)
    [Xq, Yq] = meshgrid(xq, yq);
    [Uq, Vq] = multiB(Xq, Yq, mags);                  %generating vector output (quiver)
    hold on
    plt3(X, Y, U, V, 0)                               %creating 3D line plot
    plt3(Xq, Yq, Uq, Vq, 1)
    hold off
end
