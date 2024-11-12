% Stergios Grigoriou 9564
% grigster@ece.auth.gr

%% Function for calculating low bias autocorrelation
%Inputs- 
%ts timeseries vector
%maxlag maximum lag integer
%Output-
%acf vector
%lags vector
function [acf,lags] = myautocorrelation(ts,maxlag)
    mu = mean(ts);
    s2x = var(ts,1);
    n = length(ts);
    acf = zeros(maxlag+1,1);
    lags = 0:maxlag;
    for i = lags
        acf(i+1) = sum((ts(i+1:end) - mu).*(ts(1:end-i) - mu))/(s2x*(n-i)); 
    end