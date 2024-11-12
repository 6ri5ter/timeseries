% Stergios Grigoriou 9564
% grigster@ece.auth.gr

clc
clear
close all
%% Loading data
rng(42)                    %Setting rng seed for reproducibility.
X = load("eruption2004.dat");
n = length(X);
s_index = randi(n-500,1);
X500 = X(s_index:s_index+499);
dates = {'all','500 observations'};
data = {X,X500};
acf = zeros(21,2);
pacf = acf;
lags = acf(:,1);
Q = acf(1:20,:);
Xs = Q(:,1);
for i = 1:2
    [acf(:,i),lags,Q(:,i),Xs] = stage_1_1(data{i},dates{i});
    pacf(:,i) = parcorr(data{i},NUMstd = 1.96,method = 'ols');
    figure('Name',dates{i},'NumberTitle','off')
        stem(0:20,pacf(:,i))
        hold on
        cb = 1.96/sqrt(size(data{i},1));
        yline(cb,'r--')
        yline(-cb,'r--')
        legend('PACF','95% confidence bounds')
        ylabel('PACF')
        xlabel('lag')
        title('Partial autocorrelation with 95% confidence bounds.')
        grid on
        hold off
end