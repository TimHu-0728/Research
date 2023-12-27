function [fitresult, gof] = createFit_PhiVsChi0(Chi0_data, phi_data)
%CREATEFIT(CHI0_DATA,PHI_DATA)
%  Create a fit.
%
%  Data for 'untitled fit 4' fit:
%      X Input: Chi0_data
%      Y Output: phi_data
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Dec-2023 14:40:21


%% Fit: 'untitled fit 4'.
[xData, yData] = prepareCurveData( Chi0_data, phi_data );

% Set up fittype and options.
ft = fittype( 'poly3' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 4' );
% h = plot( fitresult, xData, yData );
% legend( h, 'phi_data vs. Chi0_data', 'untitled fit 4', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'Chi0_data', 'Interpreter', 'none' );
% ylabel( 'phi_data', 'Interpreter', 'none' );
% grid on


