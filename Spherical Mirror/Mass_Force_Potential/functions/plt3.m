function plt3(X, Y, U, V, arrows)%plots
    R = sqrt(U .^ 2 + V .^ 2);
    if arrows
        qui=quiver(X, Y, U, V);
        % qui.Color = 'black'; qui.LineWidth = 1;
    else
        R = sqrt(U .^ 2 + V .^ 2);
        contourf(X, Y, log10(R), 20);
        colorbar();
    end
end