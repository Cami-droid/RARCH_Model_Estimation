% Estimation of scalar DCC(m,n) and ADCC(m,l,n) multivarate volatility model with with TARCH(p,o,q)
% or GJRGARCH(p,o,q) conditional variances

% Add 'functions' folder to the path 
    addpath('functions'); 


% USAGE:
%  [PARAMETERS] = dcc(DATA,[],M,L,N)
%  [PARAMETERS,LL,HT,VCV,SCORES,DIAGNOSTICS] = dcc(DATA,DATAASYM,M,L,N,P,O,Q,...
%                                                    GJRTYPE,METHOD,COMPOSITE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
rt=log_returns-mean(log_returns);
data=rt;
%   DATAASYM     - [OPTIONAL] K by K by T array of asymmetric covariance estimators only needed if
%                    DATA is 3-dimensional and O>0 or L>0
dataAsym=[];
%   M            - Order of symmetric innovations in DCC model
m=1;
%   L            - Order of asymmetric innovations in ADCC model
l=0;
%   N            - Order of lagged correlation in DCC model
n=1;
%   P            - [OPTIONAL] Positive, scalar integer representing the number of symmetric innovations in the
%                    univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 1.
p=1;
%   O            - [OPTIONAL] Non-negative, scalar integer representing the number of asymmetric innovations in the
%                    univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 0.
o=0;
%   Q            - [OPTIONAL] Non-negative, scalar integer representing the number of conditional covariance lags in
%                    the univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 1.
q=1;
%   GJRTYPE      - [OPTIONAL] Either 1 (TARCH/AVGARCH) or 2 (GJR-GARCH/GARCH/ARCH). Can also be a K by 1 vector
%                    containing the model type for each for each series. Default is 2.
gjrType=2;
%   METHOD       - [OPTIONAL] String, one of '3-stage' (Default) or '2-stage'.  Determines whether
%                    the model is estimated using the 3-stage estimator, or if the correlation intercepts
%                    are jointly estimated along with the dynamic parameters.
method='2-stage';
%   COMPOSITE    - [OPTIONAL] String value, either 'None' (Default), 'Diagonal' or 'Full'.  None
%                    uses standard QMLE.  'Diagonal' and 'Full' both uses composite likelihood where
%                    'Diagonal' uses all pairs of the form i,i+1 while 'Full' uses all pairs.
composite='None';
%   STARTINGVALS - [OPTIONAL] Vector of starting values to use.  See parameters and COMMENTS.
startingVals=[];
%   OPTIONS      - [OPTIONAL] Options to use in the model optimization (fmincon)

options = optimset( 'MaxIter', 30, 'MaxFunEvals', 200);
      % Limitar a 500 iteraciones
      % Limitar el n�mero m�ximo de evaluaciones de la funci�n
  
%   SPECIFICATION -[OPTIONAL] String value, either 'Scalar'(Default),'Diagonal','CP' (Common Persistence). When 'CP' is selected, l is set to zero 
specification='Diagonal';
disp('inputs values charged in the worksapace. Please run this function:');
disp('[parameters, ll ,Ht, VCV, scores, diagnostics]=rdcc(data,dataAsym,m,l,n,p,o,q,gjrType,method,composite,startingVals,options,specification)');