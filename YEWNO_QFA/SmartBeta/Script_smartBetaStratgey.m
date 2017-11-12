%% This script develops a smart beta strategy and compares the performance 
% of the portfolio built to a benchmark. The smart beta stratgey
% implemented in this code is an equal weightage stratgey, where equal
% weights are assigned to all the assets in the portfolio rather than
% weighing them according to their market capitalization and the benchmark
% considered here is SPDR
%
%% Import data from excel sheet 
prices = importfile('priceData.xlsx');
assetPrices = prices{:,2:end-1};
dates = prices(:,1);
benchMarkIndexPrices = prices{:,end};

%% Convert price data to returns
returns_portfolio = (diff(assetPrices)./assetPrices(1:end-1,:));

% Calculate mean of returns
meanReturn_portfolio = mean(returns_portfolio);

% Calculate covariance of returns
covReturns_portfolio  = cov(returns_portfolio);

%% Assign weights to each stock in the portfolio
% Equally weighted portfolio, consdiers the sum of weights assigned to all 
% assets be 100 and the weights are equal
weights = ones(1,size(returns_portfolio,2));
weights = (weights./size(returns_portfolio,2)).*100;

%% Calculate portfolio risk and return
[portfolioRisk , portfolioReturn] = portstats(meanReturn_portfolio, covReturns_portfolio , weights);

%% Return and risk of benchmark when 100 USD is invested in it
benchMark_returnsTable = (diff(benchMarkIndexPrices)./benchMarkIndexPrices(1:end-1,1));
benchMarkReturn = mean(benchMark_returnsTable)*100;
benchMarkRisk = std(benchMark_returnsTable)*100;

%% Comparing the portfolio against the benchmark
% Sharpe Ratio of the portfolio
Sharpe_Portfolio = ((portfolioReturn) - 0.0105)/portfolioRisk;

% Sharpe ratio of the benchmark
Sharpe_BenchMark = ((benchMarkReturn) - 0.0105)/benchMarkRisk;

if Sharpe_Portfolio>Sharpe_BenchMark
    disp('The portfolio outperforms the benchamrk')
else
    disp('The portfolio doesnt outperform the benchmark')
end