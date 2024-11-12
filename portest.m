% Stergios Grigoriou 9564
% grigster@ece.auth.gr

%% Function for performing portmanteau test on a time series autocorrelation
%In-
%r acf vector of a time series
%n length of the time series
%alpha significance level scalar
%plotbool plot boolean handle
%
%Out-
%Q vector of the test Q statistic for increasing lag
%Xs vector of the X^2 statistic for alpha significance level for increasing
%lag
%h test result vector
function [Q,Xs,h] = portest(r,n,alpha,plotbool)
    if nargin < 4
        plotbool = 0;
    end
    lags = 0:length(r)-1;
    nl = n-lags(2:end);
    Q = zeros(20,1);
    Xs = Q;
    h= Q;
    for i = 1:20
        Q(i) = n*(n+2)*sum((r(2:1+i).^2)./(nl(1:i)'));
        Xs(i) = chi2inv(1-alpha,i);
        h(i) = Q(i) > Xs(i);
    end
    if plotbool
        titln = ['Portmanteau test vs X^2 (', num2str((1-alpha)*100),'% confidence).'];
        figure
        plot(lags(2:end),Q)
        hold on
        plot(lags(2:end),Xs)
        legend('sample Q','X^2(k,1)')
        ylabel('Q(k)')
        xlabel('lag')
        title(titln)
    end