% Stergios Grigoriou 9564
% grigster@ece.auth.gr

close all
clear
clc
%% Loading data
rng(42)                    %Setting rng seed for reproducibility.
dates = {'1989','2000','2011'};  %Setting the requested dates and calculating how
dates_l = length(dates);   %many there are.
cdata = cell(dates_l,1); %Initializing data matrix.
s_index = [1,zeros(1,dates_l)]; % Starting index of each period
for i = 1:dates_l          %Loading data.
    dname = ['eruption',dates{i},'.dat'];
    cdata{i}= load(dname);
end
l1 = length(cdata{1});
data = [cdata{1},zeros(l1,2)];%Storing timeseries in a l1xdates_l matrix
for i = 2:dates_l
    s_index = randi(length(cdata{i})-350,1);
    data(:,i) = cdata{i}(s_index:s_index +l1-1);
end
%% Stage 1
acf = zeros(21,dates_l);
pacf = acf;
lags = acf(:,1);
Q = acf(1:20,:);
Xs = Q(:,1);
maxorder = 5; %We don't want high order models,in order to adress overffiting
keepout = 3; % How many steps ahead we will forecast
nAIC = zeros(maxorder + 1,maxorder + 1,dates_l);
%predMSE = zeros(maxorder + 1,maxorder + 1,dates_l);
indAIC = zeros(dates_l,2);
indMSE = indAIC;
%Finding 4 best models (fitting and forecast) for each year
for i = 1:dates_l
    [acf(:,i),lags,Q(:,i),Xs] = stage_1_1(data(:,i),dates{i});
    pacf(:,i) = parcorr(data(:,i),NUMstd = 1.96,method = 'ols');
    figure('Name',dates{i},'NumberTitle','off')
        stem(0:20,pacf(:,i))
        hold on
        cb = 1.96/sqrt(size(data,1));
        yline(cb,'r--')
        yline(-cb,'r--')
        legend('PACF','95% confidence bounds')
        ylabel('PACF')
        xlabel('lag')
        title('Partial autocorrelation with 95% confidence bounds.')
        grid on
        hold off
        [nAIC(:,:,i),temp,indAIC(i,:),indMSE(i,:)] = orderident(data(:,i),maxorder,keepout,1,dates{i});
        %predMSE = temp(:,:,keepout);
end
%Validate the best model with 3 step forecast 3 times in the last 9 observations.
%For 1989 after observing the plots and the data from the previous loop the
%best models seem to be AR(5),ARMA(1,2),ARMA(3,4)
%For 2000 after observing the plots and the data from the previous loop the
%best models seem to be AR(3),ARMA(1,3),ARMA(2,2)
%For 2011 after observing the plots and the data from the previous loop the
%best models seem to be AR(2),ARMA(1,2),ARMA(3,3)
order = [5,0,0;1,0,2;3,0,4];
order(:,:,2) = [3,0,0;1,0,3;2,0,2];
order(:,:,3) = [2,0,0;1,0,2;3,0,3];
predMSE = zeros(3,3);
nAIC = predMSE;
fitMSE = predMSE;
for i = 1:3
    for j = 1:3
        [predMSE(i,j),nAIC(i,j),fitMSE(i,j)] = threefoldVal(data(:,i),order(j,:,i),10);
    end
end
% Since prediction error is really close between ARMA(1,2) and ARMA(3,3) 
%for year 2011 we choose the lower order model. And we choose ARMA(1,2) and
%ARMA(1,3) for the year 1989 and 2000 respectively which achieve the lowest
% prediction error amongst the models.
order = [1,0,2;1,0,3;1,0,2];
predictions = 3;
predMSE = zeros(3,predictions);
nAIC = predMSE;
fitMSE = predMSE;
mdls = cell(3,1);
res = mdls;
for i = 1:3
    [mdls{i},res{i},~,~,fitMSE(i),predMSE(i,:),nAIC(i)] = fitARIMA(data(:,i),order(i,:),predictions,1,3,dates{i});
end
%% Portmanteau test on the residuals of the estimated models.
Qres = zeros(20,3);
figure
plot(Xs)
hold on
for i = 1:3 
    [r,la] = myautocorrelation(res{i},20);
    [Qres(:,i),~,h] = portest(r,length(res{i}),0.05);
    plot(Qres(:,i))
end
legend('X^2(k)','1989','2000','2011')
title('Portmanteau test on the residuals of the estimated models.')
xlabel('k')
ylabel('Q(k)')
