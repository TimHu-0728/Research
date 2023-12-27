function [fitresult, gof] = createFit_Ms_phi(phi_data, Ms_data)
%CREATEFIT(PHI_DATA,MS_DATA)
%  创建一个拟合。
%
%  要进行 'Ms vs. φ' 拟合的数据:
%      X 输入: phi_data
%      Y 输出: Ms_data
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 02-Nov-2023 17:51:44 自动生成


%% 拟合: 'Ms vs. φ'。
[xData, yData] = prepareCurveData( phi_data, Ms_data );

% 设置 fittype 和选项。
ft = fittype( 'poly1' );

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft );

% 绘制数据拟合图。
% figure( 'Name', 'Ms vs. φ' );
% h = plot( fitresult, xData, yData );
% legend( h, 'Ms_data vs. phi_data', 'Ms vs. φ', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % 为坐标区加标签
% xlabel( 'phi_data', 'Interpreter', 'none' );
% ylabel( 'Ms_data', 'Interpreter', 'none' );
% grid on


