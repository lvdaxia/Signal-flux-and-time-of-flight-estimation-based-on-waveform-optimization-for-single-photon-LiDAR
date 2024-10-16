function [fitresult, gof] = createFit(tra_recover_nor, Ns_recover_nor)
%CREATEFIT(TRA_RECOVER_NOR,NS_RECOVER_NOR)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : tra_recover_nor
%      Y Output: Ns_recover_nor
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 06-Apr-2024 12:01:54 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( tra_recover_nor, Ns_recover_nor );

% Set up fittype and options.
ft = fittype( 'gauss4' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 0 -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.Robust = 'LAR';
opts.StartPoint = [0.00728031541640431 218.912491690585 0.332727964809446 0.00647283338006824 219.587491690585 0.369891396507293 0.00617272709588556 218.362491690585 0.464985530935887 0.00476201690142912 220.312491690585 0.51290597245554];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'Ns_recover_nor vs. tra_recover_nor', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel tra_recover_nor
ylabel Ns_recover_nor
grid on
end


