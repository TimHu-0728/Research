function plt3(X, Y, U, V, arrows)%plots
    R       = sqrt(U .^ 2 + V .^ 2);
    if arrows                                         %choosing to plot quiver map
        qui           = quiver(X, Y, U ./ R, V ./ R); 
        qui.Color     = 'black';
        qui.LineWidth = 1.25;
    else                                              %choosing to plot contour map
        R             = sqrt(U .^ 2 + V .^ 2);        %calculating the magnitude
        contourf(X, Y, R, 20);                        %showing a range of color of 20 shades
        colorbar();                                   %displaying colorbar
    end
end
