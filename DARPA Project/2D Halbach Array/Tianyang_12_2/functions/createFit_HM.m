function [fitresult1,fitresult2,fitresult3,fitresult4,fitresult5,fitresult6,fitresult7, gof1,gof2,gof3,gof4,gof5,gof6,gof7] = createFit_HM(M20, H20, M40, H40, M60, H60, M80, H80, M100, H100, M120, H120, M140, H140)
%CREATEFIT_HM(M20, H20, M40, H40, M60, H60, M80, H80, M100, H100, M120, H120, M140, H140)
%           Create fit curves of H vs. M based on the data of Demagnetization Curve of N52 magnets 
%
%  要进行 'H vs. M ' 拟合的数据:
%      X 输入: M []
%      Y 输出: H20
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 25-Oct-2023 16:18:43 自动生成


%% 拟合: 'H vs. M (20 °C)'。
[xData20, yData20] = prepareCurveData( M20, H20 );
[xData40, yData40] = prepareCurveData( M40, H40 );
[xData60, yData60] = prepareCurveData( M60, H60 );
[xData80, yData80] = prepareCurveData( M80, H80 );
[xData100, yData100] = prepareCurveData( M100, H100 );
[xData120, yData120] = prepareCurveData( M120, H120 );
[xData140, yData140] = prepareCurveData( M140, H140 );
% 设置 fittype 和选项。
ft = 'pchipinterp';
opts = fitoptions( 'Method', 'PchipInterpolant' );
opts.ExtrapolationMethod = 'pchip';
opts.Normalize = 'on';

% 对数据进行模型拟合。
[fitresult1, gof1] = fit( xData20, yData20, ft, opts );
[fitresult2, gof2] = fit( xData40, yData40, ft, opts );
[fitresult3, gof3] = fit( xData60, yData60, ft, opts );
[fitresult4, gof4] = fit( xData80, yData80, ft, opts );
[fitresult5, gof5] = fit( xData100, yData100, ft, opts );
[fitresult6, gof6] = fit( xData120, yData120, ft, opts );
[fitresult7, gof7] = fit( xData140, yData140, ft, opts );


% 绘制数据拟合图。
figure( 'Name', 'H vs. M (Intrinsic Curve)' );
h1 = plot( fitresult1, xData20, yData20 );
h1(2).LineWidth = 1.5;
h1(2).Color     = [0 0.4470 0.7410];
h1(1).MarkerSize = 15;
h1(1).Color = [0 0.4470 0.7410];

hold on
h2 = plot( fitresult2, xData40, yData40 );
h2(2).LineWidth = 1.5;
h2(2).Color     = [0.8500 0.3250 0.0980];
h2(1).MarkerSize = 15;
h2(1).Color = [0.8500 0.3250 0.0980];

h3 = plot( fitresult3, xData60, yData60 );
h3(2).LineWidth = 1.5;
h3(2).Color     = [0.9290 0.6940 0.1250];
h3(1).MarkerSize = 15;
h3(1).Color = [0.9290 0.6940 0.1250];

h4 = plot( fitresult4, xData80, yData80 );
h4(2).LineWidth = 1.5;
h4(2).Color     = [0.4940 0.1840 0.5560];
h4(1).MarkerSize = 15;
h4(1).Color = [0.4940 0.1840 0.5560];


h5 = plot( fitresult5, xData100, yData100 );
h5(2).LineWidth = 1.5;
h5(2).Color     = [0.4660 0.6740 0.1880];
h5(1).MarkerSize = 15;
h5(1).Color = [0.4660 0.6740 0.1880];

h6 = plot( fitresult6, xData120, yData120 );
h6(2).LineWidth = 1.5;
h6(2).Color     = [0.3010 0.7450 0.9330];
h6(1).MarkerSize = 15;
h6(1).Color = [0.3010 0.7450 0.9330];

h7 = plot( fitresult7, xData140, yData140 );
h7(2).LineWidth = 1.5;
h7(2).Color     = [0.6350 0.0780 0.1840];
h7(1).MarkerSize = 15;
h7(1).Color = [0.6350 0.0780 0.1840];


% 为坐标区加标签
xlabel( 'M ($\frac{A}{m}$)', 'Interpreter', 'latex' ,'FontSize',22);
ylabel( 'H ($\frac{A}{m}$)', 'Interpreter', 'latex' ,'FontSize',22);
ylim([-10e5 0])
grid on


