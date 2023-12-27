function [fitresult, gof] = createFit_RhoVsChi0(Chi0_data, rho_data)
%CREATEFIT(CHI0_DATA,RHO_DATA)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: Chi0_data
%      Y Output: rho_data
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Dec-2023 14:12:49


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( Chi0_data, rho_data );

% Set up fittype and options.
ft = fittype( 'poly3' );
excludedPoints = excludedata( xData, yData, 'Indices', [4 11] );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'rho_data vs. Chi0_data', 'Excluded rho_data vs. Chi0_data', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'Chi0_data', 'Interpreter', 'none' );
% ylabel( 'rho_data', 'Interpreter', 'none' );
% grid on

