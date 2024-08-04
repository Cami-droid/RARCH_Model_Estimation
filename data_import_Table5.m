% Main file to run all the functions and generate the final table
clear;clc;

data = readtable('D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\data\stock_prices_28_KFT_UTX_1.csv');
% Extract data from relevant columns, Calculate log returns and Combine them in a matrix
dates = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
T=height(data)-1;
log_retrurns=zeros(1,28);
log_ret_AA=diff(log(data.AA))*100;
log_ret_AXP=diff(log(data.AXP))*100;
log_ret_BA=diff(log(data.BA))*100;
log_ret_BAC=diff(log(data.BAC))*100;
log_ret_CAT=diff(log(data.CAT))*100;
log_ret_CSCO=diff(log(data.CSCO))*100;
log_ret_CVX=diff(log(data.CVX))*100;
log_ret_DD=diff(log(data.DD))*100;
log_ret_DIS=diff(log(data.DIS))*100;
log_ret_GE=diff(log(data.GE))*100;
log_ret_HD=diff(log(data.HD))*100;
log_ret_HPQ=diff(log(data.HPQ))*100;
log_ret_IBM=diff(log(data.IBM))*100;
log_ret_INTC=diff(log(data.INTC))*100;
log_ret_JNJ=diff(log(data.JNJ))*100;
log_ret_JPM=diff(log(data.JPM))*100;
log_ret_KO=diff(log(data.KO))*100;
log_ret_MCD=diff(log(data.MCD))*100;
log_ret_MMM=diff(log(data.MMM))*100;
log_ret_MRK=diff(log(data.MRK))*100;
log_ret_MSFT=diff(log(data.MSFT))*100;
log_ret_PFE=diff(log(data.PFE))*100;
log_ret_PG=diff(log(data.PG))*100;
log_ret_T=diff(log(data.T))*100;
log_ret_TRV=diff(log(data.TRV))*100;
log_ret_VZ=diff(log(data.VZ))*100;
log_ret_WMT=diff(log(data.WMT))*100;
log_ret_XOM=diff(log(data.XOM))*100;

% Combine log returns in a matrix
log_returns = [log_ret_AA, log_ret_XOM, log_ret_BA, log_ret_BAC,log_ret_CAT log_ret_CSCO log_ret_CVX log_ret_DD log_ret_DIS log_ret_GE log_ret_HD log_ret_HPQ log_ret_IBM log_ret_INTC log_ret_JNJ log_ret_JPM log_ret_KO log_ret_MCD log_ret_MMM log_ret_MRK log_ret_MSFT log_ret_PFE log_ret_PG log_ret_T log_ret_TRV log_ret_VZ log_ret_WMT log_ret_XOM ] ;
d = size(log_returns, 2);
T = size(log_returns, 1);

% Assure that the date vector equalizes log returns' length
dates = dates(2:end); 
Task='Table5';