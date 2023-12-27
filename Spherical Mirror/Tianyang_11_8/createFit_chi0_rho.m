function [fitresult, gof] = createFit_chi0_rho(chi0, rho)
%CREATEFIT(RHO,CHI0)
%  创建一个拟合。
%
%  要进行 'ρ vs. Chi0' 拟合的数据:
%      X 输入: chi0
%      Y 输出: rho
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 24-Oct-2023 13:29:18 自动生成


%% 拟合: 'Chi0 vs. ρ'。
[xData, yData] = prepareCurveData(chi0, rho );

% 设置 fittype 和选项。
ft = fittype( 'poly1' );
excludedPoints = excludedata( xData, yData, 'Indices', 9 );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';
opts.Exclude = excludedPoints;

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% 绘制数据拟合图。
figure( 'Name', 'ρ vs. Chi0' );
h = plot( fitresult, xData, yData, excludedPoints);
h(1).MarkerSize = 10;
h(2).MarkerSize = 10;
h(2).Color=[0.8500 0.3250 0.0980];
h(3).LineWidth  = 1.5;
legend( h, '$\rho$ vs. $\chi_{0}$', 'Excluded data', 'Linear fitting curve', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% 为坐标区加标签
xlabel( '$\chi_{0}$', 'Interpreter', 'latex','FontSize',18 );
ylabel( '$\rho (\frac{kg}{m^3})$ ', 'Interpreter', 'latex','FontSize',18 );
grid on

set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
    'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [-16.8910    0.7620   16   11]); % [0 0 16 12]
set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);


