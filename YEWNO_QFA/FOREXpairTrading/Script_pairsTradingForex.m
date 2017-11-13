%% Load the data
data = importData('data.xlsx');
forex_rates = data{:,2:end};
dates = data{:,1};

% Visualzing the data of two Forex rates
figure;
plot(dates,forex_rates(:,1))
hold on
plot(dates,forex_rates(:,4))
hold on
plot(dates,forex_rates(:,10))
axis tight
title('Forex rates over time')
xlabel('Date')
ylabel('Forex Rate')
legend({'USD/JPY', 'CHFJPY', 'CADJPY'})

%% Perform pair wise analysis
numFrates = size(forex_rates,2);

reportTable = {'Forex_rate_1', 'Forex_rate_2', 'EGCI_h', 'EGCI_p', 'Corr'};

% Checking for cointegartion among each pairs of all the Forex rates
% available in the data
for i = 1:numFrates
    for j = 1:numFrates
        X = [forex_rates(:,i), forex_rates(:,j)];
        [h,p] = egcitest(X);
        ThisCorr = corr(X);
        ThisCorr = ThisCorr(1,2);
        reportTable(end+1,:) = {data.Properties.VariableNames(i+1), data.Properties.VariableNames(j+1), h, p, ThisCorr};%#ok
    end
end
warning('off')
reportTable = array2table(reportTable(2:end,:),'VariableNames', reportTable(1,:));
reportTable.EGCI_h = cell2mat(reportTable.EGCI_h);
reportTable.EGCI_p = cell2mat(reportTable.EGCI_p);
reportTable.Corr = cell2mat(reportTable.Corr);

%% Look for good pairs
% A "good" pair (h == 1)
index = find(reportTable.EGCI_h == 1 & reportTable.EGCI_p < 0.5 & reportTable.Corr > 0.8);
Pair1 = reportTable(index(2),:);
forex_rate_1 = Pair1.Forex_rate_1{1,1};
forex_rate_2 = Pair1.Forex_rate_2{1,1};

% A "bad" pair (h == 0)
Pair2 = find(reportTable.EGCI_h == 0 & reportTable.EGCI_p > 0.5 & reportTable.Corr < 0.5);

% Open a new figure window
figure
% Pair 1
N1 = find(strcmp(data.Properties.VariableNames, forex_rate_1));
N2 = find(strcmp(data.Properties.VariableNames, forex_rate_2));
subplot(2,1,1);
plot(dates, data{:,N1})
legend(data.Properties.VariableNames(N1))
title(data.Properties.VariableNames(N1))
hold on
subplot(2,1,2);
plot(dates, data{:,N2}, 'g')
legend(data.Properties.VariableNames(N2))
title(data.Properties.VariableNames(N2))
hold off

% Plotting the cointegrated series
Y = [data{:,N1},data{:,N2}];
[~,pValue,stat,cValue,reg] = egcitest(Y);
a = reg.coeff(1);
b = reg.coeff(2);
figure;
hold on
plot(dates,Y*[1;-b]-a)
legend('Cointegrated series')
title('Cointegrated series')
grid on
% Testing for stataionarity
h = adftest(Y*[1;-b]-a);
if h==1
    disp('The cointegrated series is stationary')
end

% The new cointgrated time series(spread)
spread = Y*[1;-b]-a;
stdv_spread = std(spread);
upperThreshold = 2*stdv_spread*ones(131,1);
lowerThreshold = -2*stdv_spread*ones(131,1);
figure;
plot(dates, spread);
hold on
plot(dates, upperThreshold,'b');
hold on
plot(dates, lowerThreshold,'r');
hold on
legend('spread','spreadUpperThreshold','spreadLowerThreshold')
title('Spread with its upper and lower thresholds');


%% Implementing pairs trading strategy on the testing data
% Enter and closing of position, testing it for the last 31 days in the
% data
test_spread = spread(101:131,:);

% Set a cell array which would record the positions during the trade
% Set 1 for long USDJPY and short GBPAUD and -1 for long GBPAUD and short USDJPY
position = nan(30,1);

for i = 1:length(test_spread)
    if test_spread(i) > 2*stdv_spread
        disp('long USDJPY and short GBPAUD')
        position(i) = 1;
        
    elseif test_spread(i) < -2*stdv_spread
        disp('long GBPAUD and short USDJPY')
        position(i) = -1;
    end
end

% To represent holding fill the NaN values with the previous values
position = fillmissing(position,'previous');
% if any NaN's still exists replace them with 0
for k = find(isnan(position))
  position(k) = 0;
end
%% Performance calculations
% Calculate daily market return
returns = zeros(30,1);
for i = 1: length(position)
    if position(i) == 1
        % if 100 is invested in USDJPY for every 1 in GBPAUD
        returns(i) = 100*(((Y(100+i,1)-Y(99+i,1))/Y(99+i,1))) - (((Y(100+i,2)-Y(99+i,2))/Y(99+i,2)));
    elseif position(i) == -1
        % if 100 is invested in GBPAUD for every 1 in USDJPY
        returns(i) = 100*(((Y(100+i,2)-Y(99+i,2))/Y(99+i,2))) - (((Y(100+i,1)-Y(99+i,1))/Y(99+i,1)));
    else
        returns(i) = 0;
    end
end
return_overall = sum(returns);

if return_overall > 0
    disp('winning strategy')
else
    disp('loosing stratey')
end

sharpeRatio = mean(returns)/std(returns);
disp('The sharpeRatio is')
disp(sharpeRatio)
