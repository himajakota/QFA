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
index = find(reportTable.EGCI_h == 1 & reportTable.EGCI_p < 0.5);
Pair1 = reportTable(index(2),:);
forex_rate_1 = Pair1.Forex_rate_1{1,1};
forex_rate_2 = Pair1.Forex_rate_2{1,1};
 
% A "bad" pair (h == 0)
Pair2 = find(reportTable.EGCI_h == 0 & reportTable.EGCI_p > 0.5 & reportTable.Corr < 0.5);

% open a new figure window
figure
% Pair 1
N1 = find(strcmp(data.Properties.VariableNames, forex_rate_1));
N2 = find(strcmp(data.Properties.VariableNames, forex_rate_2));
subplot(2,1,1);
plot(data.Date, data{:,N1})
legend(data.Properties.VariableNames(N1))
title(data.Properties.VariableNames(N1))
hold on
subplot(2,1,2); 
plot(data.Date, data{:,N2}, 'g')
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
% The new cointgrated time series 
ts = Y*[1;-b]-a;
if h==1
    disp('The cointegrated series is stationary')
end
