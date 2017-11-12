%% This script predicts GDP(economic indicator) using a time series ARIMA model

%% Import data
gdp_data = importData('GDP.csv');
date = gdp_data{:,1};
data = gdp_data{:,2};

% Visulaization of the time series over the time period
figure;
plot(date,data);
hold on;
title('GDP');
xlabel('time');
ylabel('GDP');
legend('GDP')
hold off;

% Test for stationarity of the data
h = adftest(data);
if h == false
disp('Fail to reject the null hypothesis implies that GDP is not trend stationary.')
end
%% The time series has an upward trend. Check the stationairty of the series by taking its first difference
% Create a differencing lag operator polynomial object, and then use it to filter the observed series.
D1 = LagOp({1,-1},'Lags',[0,1]);
diffData = filter(D1,data);
N = length(data);

% Visualizing the differenced data
figure
plot(2:N,diffData)
xlim([0,N])
legend('First differenced GDP')
title('First Differenced GDP Series')
% Test for stationarity of the differenced data
h1 = adftest(diffData); % as the null hypothesis was rejected move ahead to find the parameters of ARIMA model
if h1 == true
disp('Rejection of the null implies tha GDP is trend stationary.')
end
%% Plot ACF and PACF
figure
subplot(2,1,1)
autocorr(diffData);
subplot(2,1,2)
parcorr(diffData);

%% Creating the ARIMA model
% Estimating the best parameter values for ARIMA model,with the help of information
% criteria,use a combination various lags for both ar and ma, a total of 16
% models and select the model with the least AIC value 
maxLags = 4;
AIC_ofAllModels = zeros(maxLags,maxLags);
parfor i = 1:maxLags
    for j = 1:maxLags
        mdl = arima(i,1,j);
        [~, ~, logL, info] = estimate(mdl, data, 'display', 'off');
        AIC_ofAllModels(i, j) = aicbic(logL, length(info.X));
    end
end

[arLag, maLag] = find(AIC_ofAllModels==min(min(AIC_ofAllModels)));

% Build the model with the optimized ar and ma lags
mdl = arima(arLag, 1, maLag);
EstMdl = estimate(mdl,data);

% Diagnostic checking of the residual, it should a white noise
res = infer(EstMdl,data);

% Plot the residual
figure
subplot(2,2,1)
plot(res./sqrt(EstMdl.Variance))
title('Standardized Residuals')
subplot(2,2,2)
qqplot(res)
subplot(2,2,3)
autocorr(res,8)
subplot(2,2,4)
parcorr(res,8)

%% Generate forecast for the next 5 years
[dataPredicted,MSE] = forecast(EstMdl,20,'Y0',data);
upperBound = dataPredicted + 1.96*sqrt(MSE);
lowerBound = dataPredicted - 1.96*sqrt(MSE);
T = length(date);

figure
plot(data);
hold on
plot(284:303,dataPredicted,'r','LineWidth',2);
plot(284:303,upperBound,'b--','LineWidth',1.5);
plot(284:303,lowerBound,'b--','LineWidth',1.5);
futureDates = [date; date(T) + cumsum(diff(date(T-20:T)))];
h = gca;
h.XTick = 1:10:(T+16);
h.XTickLabel = year(datestr(futureDates(1:10:end),21));
legend('GDP');
title('GDP forecast for 5 years')
hold off
