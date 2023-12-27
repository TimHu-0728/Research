function [fitresult, gof] = createFit_TH(H_demag, T)
%CREATEFIT(H_DEMAG,T)
%  Create a fit curve for T vs H_demag
%
%  要进行 'T vs H_demag' 拟合的数据:
%      X input: H_demag
%      Y output: T
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 25-Oct-2023 17:59:14 自动生成


%% 拟合: 'T vs H_demag'。
[xData, yData] = prepareCurveData( H_demag, T );

% 设置 fittype 和选项。
ft = 'pchipinterp';
opts = fitoptions( 'Method', 'PchipInterpolant' );
opts.ExtrapolationMethod = 'pchip';
opts.Normalize = 'on';

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% 绘制数据拟合图。
figure( 'Name', 'T vs H_demag' );
h = plot( fitresult, xData, yData );
h(2).LineWidth  = 1.5;
h(1).MarkerSize = 15;
hold on

% 为坐标区加标签
xlabel( 'H ($\frac{A}{m}$)', 'Interpreter', 'latex' ,'FontSize',22);
ylabel( 'T (°C)', 'Interpreter', 'none','FontSize',22 );
grid on


