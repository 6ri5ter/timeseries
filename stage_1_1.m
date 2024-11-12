function [acf,lags,Q,Xs] = stage_1_1(X,name)
    mu = mean(X);
    s2x = var(X,1);
    n = length(X);
    %Plotting time series
    titln = ['Wating time. (mu = ',num2str(mu),')'];
    fname = ['Plots for year ',name,'.'];
    figure('Name',fname,'NumberTitle','off')
    tiledlayout(3,1)
    nexttile
    plot(X)
    title(titln)
    titln = 'Sample autocorrelation with 95% confidence bounds.';
    %Calculating and plotting low bias acf with 95% confidence bounds for 20 lags + r(0) = 1.
    acf = zeros(21,1);
    lags = 0:20;
    for i = lags
        acf(i+1) = sum((X(i+1:end) - mu).*(X(1:end-i) - mu))/(s2x*(n-i)); 
    end
    nexttile
    stem(0:20,acf)
    hold on
    cb = 1.96/sqrt(n);
    yline(cb,'r--')
    yline(-cb,'r--')
    legend('ACF','95% confidence bounds')
    ylabel('ACF')
    xlabel('lag')
    title(titln)
    grid on
    hold off
    %Portmanteau test (Ljung-Box variation)
    nl = n-lags(2:end);
    Q = zeros(20,1);
    Xs = Q;
    for i = 1:20
        Q(i) = n*(n+2)*sum((acf(2:1+i).^2)./(nl(1:i)'));
        Xs(i) = chi2inv(0.95,i);
    end
    titln = 'Portmanteau test vs X^2 (95% confidence).';
    nexttile
    plot(lags(2:end),Q)
    hold on
    plot(lags(2:end),Xs)
    legend('sample Q','X^2(k,1)')
    ylabel('Q(k)')
    xlabel('lag')
    title(titln)