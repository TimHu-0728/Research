function  Plotting_result(H_max,data)
%PLOTTING_RESULT 
% Plot the Demagnetization curves (intrinsic curve and normal curve) of N52
% magnet in SI unit. Then plot the H_max found around the Halbach Array on
% the same figure

%%
% Temperature: 20°C
H20  = data(1,:);                          % [A/m]
B20  = data(2,:);                          % [T]

% Temperature: 40°C
H40 = data(3,:);
B40 = data(4,:);

% Temperature: 60°C
H60 = data(5,:);
B60 = data(6,:);

% Temperature: 80°C
H80 = data(7,:);
B80 = data(8,:);

% Temperature: 100°C
H100 = data(9,:);
B100 = data(10,:);

% Temperature: 120°C
H120 = data(11,:);
B120 = data(12,:);

% Temperature: 140°C
H140 = data(13,:);
B140 = data(14,:);

% Plot H_max on the graph
B    = linspace(0,1.5,100);
H    = H_max * ones(1,100);

figure
plot(H20,B20,H40,B40,H60,B60,H80,B80,H100,B100,H120,B120,H140,B140,'LineWidth',2)
hold on
plot(H,B,'k--',LineWidth=2)
title('Demagnetization curve for N52 Magnets','FontSize',24)
xlabel('H (A/m)','FontSize',18)
ylabel('B (T)','FontSize',18)
leg=legend('20 °C','40 °C','60 °C','80 °C','100 °C','120 °C','140 °C','$H_{max}$',Location='eastoutside',fontsize=18);
set(leg,'Interpreter','latex')
grid on

