function [fitresult, gof] = createFit_phi_chi0(Chi0_data, phi_data)
%CREATEFIT(CHI0_DATA,PHI_DATA)
%  创建一个拟合。
%
%  要进行 'φ vs chi0' 拟合的数据:
%      X 输入: Chi0_data
%      Y 输出: phi_data
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 02-Nov-2023 17:52:55 自动生成


%% 拟合: 'φ vs chi0'。
[xData, yData] = prepareCurveData( Chi0_data, phi_data );

% 设置 fittype 和选项。
ft = fittype( 'poly1' );
excludedPoints = excludedata( xData, yData, 'Indices', 9 );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Exclude = excludedPoints;

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% % 绘制数据拟合图。
% figure( 'Name', 'φ vs chi0' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'phi_data vs. Chi0_data', '已排除 phi_data vs. Chi0_data', 'φ vs chi0', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % 为坐标区加标签
% xlabel( 'Chi0_data', 'Interpreter', 'none' );
% ylabel( 'phi_data', 'Interpreter', 'none' );
% grid on


