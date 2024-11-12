% Stergios Grigoriou 9564
% grigster@ece.auth.gr

%% Function for helping in finding the optimal ARMA model order
%X time series vector
%maxorder maximum p and q
%keepout how many steps ahead to predict
%plotbool plotting boolean handle
%name year for the plot
%OUT
%nAIC normalised AIC matrix for each fitted model
%predMSE -//-
%indAIC/MSE index of the minimum for each case for automation
function [nAIC,predMSE,indAIC,indMSE] = orderident(X,maxorder,keepout,plotbool,name)
    if nargin < 4
        plotbool = 0;
    end
    nAIC = zeros(maxorder + 1);
    predMSE = zeros(maxorder + 1,maxorder+1,keepout);
    minAIC = inf;
    minpredMSE = inf;
    indAIC = [0,0];
    for i = 0:maxorder
        for j = 0:maxorder
            if~(j || i)
                continue
            end
            [~,~,~,~,~,predMSE(i+1,j+1,:),nAIC(i+1,j+1)] = fitARIMA(X,[i,0,j],keepout);
            if nAIC(i+1,j+1) < minAIC
                minAIC = nAIC(i+1,j+1);
                indAIC = [i,j];
            end
            if predMSE(i+1,j+1,keepout) < minpredMSE
                minpredMSE = predMSE(i+1,j+1,keepout);
                indMSE = [i,j];
            end
        end
    end
    if plotbool
        figure('Name',name)
        tiledlayout(1,2)
        nexttile
        plot(0:maxorder,nAIC(1:maxorder+1,:))
        labels = cell(1,maxorder);
        for i = 1:maxorder+1
            labels{i} = ['p = ',num2str(i-1)];
        end
        legend(labels)
        title('AIC vs order q of MA and p of AR')
        xlabel('q')
        ylabel('nAIC')
        nexttile
        plot(0:maxorder,predMSE(1:maxorder+1,:,keepout))
        labels = cell(1,maxorder);
        for i = 1:maxorder+1
            labels{i} = ['p = ',num2str(i-1)];
        end
        legend(labels)
        title('MSE of forecast 3 steps ahead vs order q of MA and p of AR')
        xlabel('q')
        ylabel('MSE')
    end
       