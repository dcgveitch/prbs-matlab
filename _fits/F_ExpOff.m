function [fitresult, gof] = createFit(X, Y)
%CREATEFIT(X,YTEST)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : X
%      Y Output: Ytest
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 15-May-2013 09:37:48


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( X, Y );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [0 -Inf -Inf];
opts.StartPoint = [5000 1.5 200];
opts.Upper = [100000 Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'Ytest vs. X', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel( 'X' );
% ylabel( 'Ytest' );
% grid on


