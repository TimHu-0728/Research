function [fitresult, gof] = createFit_Ms_rho(rho, Ms)
%CREATEFIT(RHO,MS)
%  创建一个拟合。
%
%  要进行 'Ms vs. ρ' 拟合的数据:
%      X 输入: rho
%      Y 输出: Ms
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 24-Oct-2023 13:23:26 自动生成


%% 拟合: 'Ms vs. ρ'。
[xData, yData] = prepareCurveData( rho, Ms );

% 设置 fittype 和选项。
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% 绘制数据拟合图。
% figure( 'Name', 'Ms vs. ρ' );
% h = plot( fitresult, xData, yData );
% h(1).MarkerSize = 10;
% h(2).LineWidth  = 1.5;
% legend( h, 'Ms vs. $\rho$', 'Liear Fitting Curve', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% % 为坐标区加标签
% xlabel( '$\rho (\frac{kg}{m^3})$', 'Interpreter', 'latex','FontSize',18 );
% ylabel( 'Ms', 'Interpreter', 'none','FontSize',18 );
% grid on


