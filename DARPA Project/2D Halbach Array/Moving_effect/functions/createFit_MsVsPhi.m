function [fitresult, gof] = createFit_MsVsPhi(phi_data, Ms_data)
%CREATEFIT(PHI_DATA,MS_DATA)
%  Create a fit.
%
%  Data for 'untitled fit 3' fit:
%      X Input: phi_data
%      Y Output: Ms_data
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Dec-2023 14:36:58


%% Fit: 'untitled fit 3'.
[xData, yData] = prepareCurveData( phi_data, Ms_data );

% Set up fittype and options.
ft = fittype( 'poly3' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 3' );
% h = plot( fitresult, xData, yData );
% legend( h, 'Ms_data vs. phi_data', 'untitled fit 3', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'phi_data', 'Interpreter', 'none' );
% ylabel( 'Ms_data', 'Interpreter', 'none' );
% grid on


