function plot_HBArray(X,Y,Z,Bx,By,Bz,r_mag,orient_mag,dimensions,plot_mag)
%PLOT_HBARRAY 此处显示有关此函数的摘要
%   此处显示详细说明

%%
figure
hold on
axis equal
grid on

% Plot coordinate axes arrows (x, y, z)
quiver3(0,0,0,1,0,0,0.8*dimensions(3),'r', 'LineWidth', 1,'MaxHeadSize',1.5)
quiver3(0,0,0,0,1,0,0.8*dimensions(3),'g', 'LineWidth', 1,'MaxHeadSize',1.5)
quiver3(0,0,0,0,0,1,0.8*dimensions(3),'b', 'LineWidth', 1,'MaxHeadSize',1.5)

% Plot B field
quiver3(X,Y,Z,Bx,By,Bz,'LineWidth',0.75,'MaxHeadSize',2,'Color',[0.2 0.6 0.2]);

xlabel('X (m)','FontSize',20)
ylabel('Y (m)','FontSize',20)
zlabel('Z (m)','FontSize',20)
title('B field of Halbach Array (Tesla)','FontSize',22)

% Plot magnets
n_mag = size(r_mag,1);
for i = 1:n_mag 
    plot_magnet(r_mag(i,:),orient_mag(i,:),dimensions,plot_mag)
end
legend('x','y','z','B(x,y,z)','Location','bestoutside',fontsize=18)
hold off
end

