% Stergios Grigoriou 9564
% grigster@ece.auth.gr

%% Function for validating the best models
%X time series vector
%order in the form of [p,d,q] for the best models
%keepout how many steps ahead to forecast
%OUT
%metrics of the fitted model averaged 
function [predMSE,nAIC,fitMSE] = threefoldVal(X,order,keepout)
    n = length(X);
    mu = mean(X);
    Xmu = X - mu;
    predMSE = 0;
    fitMSE = 0;
    nAIC = 0;
    for j = keepout:-1:1
        Xtrain = Xmu(1:end-keepout*j);
        %Holding some validation data for reproducible results.~(90-10)
        Xtest = Xmu(end-keepout*j + 1:end);
        mdl = arima(order(1),order(2),order(3));
        mdlest = estimate(mdl,Xtrain,'Display','off');
        res = infer(mdlest,(Xtrain));
        fitMSE = fitMSE + var(res,1);
        Xest = forecast(mdlest,keepout,Xtrain);
        predMSE = predMSE + var(Xtest(1:keepout) - Xest,1);%they have both zero mean
        nAIC = nAIC + log(mdlest.Variance) + 2*sum(order)/n;
    end
    fitMSE = fitMSE/keepout;
    predMSE = predMSE/keepout;
    nAIC = nAIC/keepout;