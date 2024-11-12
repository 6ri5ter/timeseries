% Stergios Grigoriou 9564
% grigster@ece.auth.gr

%% Function for estimating ARIMA models
% Inputs-
% X time series vector
% order vector in the form [p,d,q]
% keepout max step ahead for prediction must be integer
% plotbool true for plots false for no plots
% testl length of test set
% year charachter for the plot
% Outputs-
% mdlest the estimated arima model
% res the residuals of the fitted arima model
% Xest the estimated future values of the timeseries X
% fitMSE fitting mean squared error of the model
% predMSE prediction mean squared error of the predictions 
% nAIC normalised AIC for the model

function [mdlest,res,Xest,Xtest,fitMSE,predMSE,nAIC] = fitARIMA(X,order,keepout,plotbool,testl,name)
    if nargin < 4
        plotbool = 0;
        testl = 30;
    elseif nargin < 5
        testl = 30;
    end
    mu = mean(X);
    Xmu = X - mu;
    n = length(X);
    Xtrain = Xmu(1:end-testl);
    %Holding some validation data for reproducible results.~(90-10)
    Xtest = Xmu(end-testl+1:end);
    mdl = arima(order(1),order(2),order(3));
    mdlest = estimate(mdl,Xtrain,'Display','off');
    res = infer(mdlest,(Xtrain));
    fitMSE = var(res,1);
    Xhat = Xtrain - res + mu;
    predMSE = zeros(keepout,1);
    for i = 1:keepout
        %We calculate the forcast for upt to keepout steps ahead and the
        %corresponding mse but we only keep the max keepout estimates for
        %reference
        Xest = forecast(mdlest,i,Xtrain);
        predMSE(i) = var(Xtest(1:i) - Xest,1);%they have both zero mean
        if i == 1
            predMSE(1) = (Xtest(1) - Xest)^2;
        end
    end
    Xest = [Xhat(end);Xest + mu];
    nAIC = log(mdlest.Variance) + 2*sum(order)/n;
    Xtest = Xtest +mu;
    if plotbool
        figure
        plot(X(1:end-testl+keepout))
        hold on
        plot(Xhat)
        plot(n-testl:n-testl+keepout,Xest)
        titln = ['Observed waiting time vs ARMA estimation for year ',name,'.'];
        title(titln)
        plot(n-testl:n-testl+keepout,1.96*sqrt(mdlest.Variance) + Xest,'g--')
        plot(n-testl:n-testl+keepout,-1.96*sqrt(mdlest.Variance) + Xest,'g--')
        labels = {'Waiting time',['ARMA(',num2str(order(1)),',',num2str(order(3)),')'],'Forecasted','Forecast bounds'};
        legend(labels)
    end