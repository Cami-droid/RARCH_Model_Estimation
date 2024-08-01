% Main file to run all the functions and generate the final table
clear;clc;

data = readtable('D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\data\stock_prices_28_KFT_UTX_1.csv');
height(data);
% Extract data from relevant columns, Calculate log returns and Combine them in a matrix
dates = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
T=height(data)-1;
log_retrurns=zeros(1,10);
log_ret_BAC=diff(log(data.BAC))*100;%%%%%%%%%%%%%
log_ret_JPM=diff(log(data.JPM))*100;%%%%%%%%%%%%%%
log_ret_IBM=diff(log(data.IBM))*100;%%%%%%%%%%%%%
log_ret_MSFT=diff(log(data.MSFT))*100;%%%%%%%%%%%
log_ret_XOM=diff(log(data.XOM))*100;%%%%%%%%%%%%
log_ret_AA=diff(log(data.AA))*100;%%%%%%%%
log_ret_AXP=diff(log(data.AXP))*100;%%%%%%%%%, 
log_ret_DD=diff(log(data.DD))*100;%%%%%%%%%%%%
log_ret_GE=diff(log(data.GE))*100;%%%%%%%%%%%
log_ret_KO=diff(log(data.KO))*100;%%%%%%%%%%%%%%%

% Combine log returns in a matrix

log_returns = [log_ret_BAC,log_ret_JPM,log_ret_IBM,log_ret_MSFT,log_ret_XOM,log_ret_AA,log_ret_AXP, log_ret_DD ,log_ret_GE , log_ret_KO] ;
d = size(log_returns, 2);
T = size(log_returns, 1);

% Assure that the date vector equalizes log returns' length
dates = dates(2:end); 
Task='Table4';