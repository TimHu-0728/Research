function plot_magnetizations(X_layer,Mc,Ms,M)
    figure
    subplot(2,2,1)
    plot(X_layer, Mc(:,1), X_layer, Ms(1)*ones(size(X_layer)), X_layer, M(:,1),LineWidth=1.5)
    xlabel('X (m)','FontSize',16)
    ylabel('$M (\frac{A}{m})$','FontSize',16,'Interpreter','latex')
    ylim([2000 10000])
    legend('Mc (Critical Magnetization)','Ms (Saturation Magnetization)','M (Ferrofluid Magnetization)','Location','east',fontsize=12)
    title('$\chi_0 = 0.1$','FontSize',20,'Interpreter','latex')
    grid on
    
    subplot(2,2,2)
    plot(X_layer, Mc(:,2), X_layer, Ms(2)*ones(size(X_layer)), X_layer, M(:,2),LineWidth=1.5)
    xlabel('X (m)','FontSize',16)
    ylabel('$M (\frac{A}{m})$','FontSize',16,'Interpreter','latex')
    ylim([2000 10000])
    title('$\chi_0 = 1$','FontSize',20,'Interpreter','latex')
    grid on
    
    subplot(2,2,3)
    plot(X_layer, Mc(:,3), X_layer, Ms(3)*ones(size(X_layer)), X_layer, M(:,3),LineWidth=1.5)
    xlabel('X (m)','FontSize',16)
    ylabel('$M (\frac{A}{m})$','FontSize',16,'Interpreter','latex')
    ylim([2000 10000])
    title('$\chi_0 = 1.2$','FontSize',20,'Interpreter','latex')
    grid on
    
    subplot(2,2,4)
    plot(X_layer, Mc(:,4), X_layer, Ms(4)*ones(size(X_layer)), X_layer, M(:,4),LineWidth=1.5)
    xlabel('X (m)','FontSize',16)
    ylabel('$M (\frac{A}{m})$','FontSize',16,'Interpreter','latex')
    ylim([2000 10000])
    title('$\chi_0 = 1.5$','FontSize',20,'Interpreter','latex')
    grid on
end

