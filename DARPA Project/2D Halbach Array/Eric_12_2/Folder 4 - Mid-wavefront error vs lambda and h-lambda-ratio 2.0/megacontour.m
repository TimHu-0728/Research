function megacontour(x, y, Z, c0, c1, color)
%MEGACONTOUR This function plots a contour inequality by plotting a solid line denoting the equality, and a dashed line denoting the side of the equality you
%want to be on.

%Inputs are the x, y (note these two should be 1D), and 2D Z values for the contour plot, the equality Z = c0, and a c1 on the right side of Z = c0 to show
%where the inequality should be. color is the color in a string
    cc0 = contourc(x, y, Z,[c0, c0]);
    cc1 = contourc(x, y, Z,[c1, c1]);
    hold on
    plot(cc0(1,2:end), cc0(2,2:end), 'Color', color, 'LineStyle', '-', 'LineWidth', 3)
    plot(cc1(1,2:end), cc1(2,2:end), 'Color', color, 'LineStyle', '-.', 'LineWidth', 3)
end

