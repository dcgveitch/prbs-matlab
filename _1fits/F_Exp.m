function [fitresult, gof] = createFit(X, Ytest)
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

%  Auto-generated by MATLAB on 15-May-2013 09:09:08


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( X, Ytest );

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf];
opts.StartPoint = [1211.90179923893 -1.65060766975366];
opts.Upper = [Inf Inf];

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


