% Main file to run all the functions and generate the final table
clearvars -except table_number;

;clc;

data = readtable('D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\data\stock_prices_28_KFT_UTX_1.csv');

% Extract data from relevant columns
dates = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
AA = data.AA; 
XOM = data.XOM;

% Calculate log returns
log_returns_AA = diff(log(AA)) * 100;
log_returns_XOM = diff(log(XOM)) * 100;

% Combine log returns in a matrix
log_returns = [log_returns_AA, log_returns_XOM];

d = size(log_returns, 2);
T = size(log_returns, 1);

% Assure that the date vector equalizes log returns' length
dates = dates(2:end);
Task='Table2';