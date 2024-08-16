function [parameters, ll ,Ht, VCV, scores]=rdcc(data,m,n,p,q,method,composite,startingVals,options,specification)
% Estimation of scalar DCC(m,n) and ADCC(m,l,n) multivarate volatility model with with TARCH(p,o,q)
% or GJRGARCH(p,o,q) conditional variances
%
% USAGE:
%  [PARAMETERS] = dcc(DATA,[],M,L,N)
%  [PARAMETERS,LL,HT,VCV,SCORES,DIAGNOSTICS] = dcc(DATA,DATAASYM,M,L,N,P,O,Q,...
%                                                    GJRTYPE,METHOD,COMPOSITE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   M            - Order of symmetric innovations in DCC model
%   N            - Order of lagged correlation in DCC model
%   P            - [OPTIONAL] Positive, scalar integer representing the number of symmetric innovations in the
%                    univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 1.
%   Q            - [OPTIONAL] Non-negative, scalar integer representing the number of conditional covariance lags in
%                    the univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 1.
%   METHOD       - [OPTIONAL] String, one of '3-stage' (Default) or '2-stage'.  Determines whether
%                    the model is estimated using the 3-stage estimator, or if the correlation intercepts
%                    are jointly estimated along with the dynamic parameters.
%   COMPOSITE    - [OPTIONAL] String value, either 'None' (Default), 'Diagonal' or 'Full'.  None
%                    uses standard QMLE.  'Diagonal' and 'Full' both uses composite likelihood where
%                    'Diagonal' uses all pairs of the form i,i+1 while 'Full' uses all pairs.
%   STARTINGVALS - [OPTIONAL] Vector of starting values to use.  See parameters and COMMENTS.
%   OPTIONS      - [OPTIONAL] Options to use in the model optimization (fmincon)
%   SPECIFICATION -'Scalar','Diagonal','CP'
%
% OUTPUTS:
%   PARAMETERS   - Estimated parameters.  Output depends on METHOD.
%                    3-stage: [VOL(1) ... VOL(K) corr_vech(R)' vech(N)' alpha gamma beta]
%                    2-stage: [VOL(1) ... VOL(K) corr_vech(R)' alpha gamma beta]
%                    where VOL(j) is a (1+P(i)+O(i)+Q(i)) vector containing the parameters from
%                    volatility model i.
%   LL           - The log likelihood at the optimum
%   HT           - A [K K T] dimension matrix of conditional covariances
%   VCV          - A numParams^2 square matrix of robust parameter covariances (A^(-1)*B*A^(-1)/T)
%   SCORES       - A T by numParams matrix of individual scores
%   DIAGNOSTICS  - A structure containing further outputs.
%
% COMMENTS:
%   The dynamics of a the correlations in a DCC model are:
%     3-stage:
%     Q(t) = R*(1-sum(a)-sum(b))-sum(g)*N + a(1)*e(t-1)'*e(t-1) + ... + a(m)*e(t-m)'*e(t-m)
%     + g(1)*v(t-1)'*v(t-1) + ... + g(l)*v(t-l)*v(t-l) + b(1)*Q(t-1) + ... + b(n)*Q(t-1)
%
%     2-stage
%     Q(t) = R.*scale + a(1)*e(t-1)'*e(t-1) + ... + a(m)*e(t-m)'*e(t-m)
%     + g(1)*v(t-1)'*v(t-1) + ... + g(l)*v(t-l)*v(t-l) + b(1)*Q(t-1) + ... + b(n)*Q(t-1)
%
%   where v(t,:) = e(t,;).*(e(t,:)<0) and s = sqrt((1-sum(a)-sum(b)-[]*sum(g))) and scale = s*s'
%
%
% EXAMPLES:
%   % DCC(1,1)
%   parameters = dcc(data,[],1,0,1)
%   % ADCC(1,1)
%   parameters = dcc(data,[],1,1,1)
%   % ADCC(1,1), 2-stage
%   parameters = dcc(data,[],1,1,1,[],[],[],[],'2-stage')
%
% See also CCC_MVGARCH, BEKK, RARCH, SCALAR_VT_VECH, MATRIX_GARCH, TARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/13/2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 3
        p = [];
        q = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
        specification = [];
    case 4
        q = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
        specification = [];
    case 5
        method=[];
        composite = [];
        startingVals = [];
        options = [];
        specification = [];
    case 6
        composite=[];
        startingVals = [];
        options = [];
        specification = [];
    case 7
        startingVals=[];
        options = [];
        specification = [];
    case 8
        options = [];
        specification = [];
    case 9
        specification = [];
    case 10
        % perfect
    otherwise
        error('3 to 10 inputs required.')
end

if ismatrix(data)
    [T,k] = size(data);
    data2d = data;
   data = zeros(k,k,T);
    for t=1:T
        data(:,:,t) = data2d(t,:)'*data2d(t,:);
    end
elseif ndims(data)==3
    [k,~,T]=size(data);
    data2d = zeros(k,T);
else
    error('DATA must be either a K by T matrix of a K by K by T 3-dimensional array.')
end

if ~isscalar(m) || floor(m)~=m || m<1
    error('M must be a positive integer.')
end
if isempty(n)
    n = 0;
end
if ~isscalar(n) || floor(n)~=n || n<0
    error('N must be a non-negative integer.')
end

if isempty(p)
    p = 1;
end
if isscalar(p)
    p = p * ones(k,1);
end
if isempty(q)
    q = 1;
end
if isscalar(q)
    q = q * ones(k,1);
end
if any(floor(p)~=p) || any(p<1) || numel(p)~=k
    error('All elements in P must be positive, and P must be either scalar or K by 1')
end
if any(floor(q)~=q) || any(q<0) || numel(q)~=k
    error('All elements in Q must be non-negative, and Q must be either scalar or K by 1.')
end
gjrType = ones(k,1)*2;

if isempty(method)
    method = '3-stage';
end
method = lower(method);
if ~ismember(method,{'3-stage','2-stage'})
    error('TYPE must be either ''3-stage'' or ''2-stage''.')
end
if strcmpi(method,'3-stage')
    stage = 3;
else
    stage = 2;
end


if isempty(composite)
    composite = 'none';
end
composite = lower(composite);
if ~ismember(composite,{'none','diagonal','full'})
    error('COMPOSITE must be one of ''None'', ''Diagonal'' or ''Full''.')
end
if strcmpi(composite,'none')
    composite = 0;
elseif strcmpi(composite,'diagonal')
    composite = 1;
else
    composite = 2;
end
if stage == 2 && composite==1
    warning('oxfordMFE:incorrectOption','When TYPE is ''2-stage'', COMPOSITE must be either ''None'' or ''Full''.')
    composite = 2;
end

if ~isempty(startingVals)
    count = k + sum(p)+ sum(q); % idxM
    count = count + m + n; %idxM+%idxD
    if stage==2 
        count = count + k*(k-1)/2;  %idxM+%idxD+%idxS 2-stage
    elseif  stage==3
        count = count + k*(k+1)/2 + k*(k-1)/2; %idxM+%idxD+%idxS 3-stage
    end
    if length(startingVals)~=count
        error('STARTINGVALS does not contain the correct number of parameters.')
    end
    count = k+sum(p)+sum(q);
    tarchStartingVals = startingVals(1:count); % starting_thetaM
    offset = count;
    % Count for intercepts
    if stage==2 
        count = k*(k-1)/2;
    elseif  stage==3
        count = k*(k+1)/2 + k*(k-1)/2;
    end
    offset = offset + count;
    count = m+n;
    dccStartingVals = startingVals(offset + (1:count)); % starting_thetaD
else
    tarchStartingVals = []; % starting_thetaM
    dccStartingVals = []; %starting_thetaM
end


if isempty(options)
    options = optimset('fmincon');
    options.Display = 'iter';
    options.Diagnostics = 'on';
    options.Algorithm = 'interior-point';
end
try
    optimset(options);
catch ME
    error('OPTIONS does not appear to be a valid options structure.')
end
%added
if isempty(specification)
specification='Scalar';
end
%added-end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Univariate volatility models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o=zeros(size(p));
[H,univariate] = dcc_fit_variance(data2d,p,o,q,gjrType,tarchStartingVals); %starting_thetaM
stdData = data;
for t=1:T
    h = sqrt(H(t,:));
    hh = h'*h;
    stdData(:,:,t) = stdData(:,:,t)./hh;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Back casts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = .06*.94.^(0:sqrt(T));
w = w'/sum(w);
backCast = zeros(k);
for i=1:length(w)
    backCast = backCast + w(i)*stdData(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = mean(stdData,3); % H_bar or unconditional covariance matrix
r = sqrt(diag(R));
R = R ./ (r*r'); % Ct unconditional correlation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(dccStartingVals)
    startingVals = dccStartingVals;
else
    a = [.01 .03 .05 .1];
    theta = [.99 .97 .95];
        [a,theta] = ndgrid(a,theta);
        parameters = unique([a(:) theta(:)-a(:)],'rows');
        isAsym = 0;
    end
    minLL= inf;
    isJoint = false; % doesn't include volatility parameters
    isInference = false;
    for i=1:size(parameters,1)
        ll = dcc_likelihood(parameters(i,:),stdData,[],1,isAsym,1,R,[],backCast,[],3,composite,isJoint,isInference);
        if ll<minLL
            startingVals = parameters(i,:);
            minLL = ll;
        end
    end
    a = startingVals(1);
    b = startingVals(2);
    startingVals = [a*ones(1,m)/m b*ones(1,n)/n];
        %added
    switch specification
       case 'Diagonal'
            startingVals = [startingVals(1:m)*ones(1,k)/k , startingVals(1+m:m+n)*ones(1,k)/k];
       case 'CP'
            startingVals = [startingVals(1:m)*ones(1,k)/k, (a+b)/k];
    end
    
    

LB = zeros(length(startingVals),1);
UB = ones(length(startingVals),1);
switch  specification
case 'Scalar'
A=horzcat(kron(eye(1),ones(1,m)),kron(eye(1),ones(1,n)));
b = 0.99998;

case 'Diagonal'
A=horzcat(kron(eye(k),ones(1,m)),kron(eye(k),ones(1,n)));
b = repmat(0.99998, k, 1);

case 'CP'
A=horzcat(kron(eye(k),ones(1,m)),ones(k,1));
b = repmat(0.99998, k, 1);
end


%if (startingVals*A'-b) >= 0
 %   error('STARTINGVALS for DCC parameters are not comparible with a positive definite intercept.')
%end

isJoint = false; % doesn't include volatility parameters
isInference = false;
parameters = fmincon(@rdcc_likelihood,startingVals,A,b,[],[],LB,UB,[],options,stdData,m,n,R,backCast,3,composite,isJoint,isInference,[],specification);

switch specification
    case 'Scalar'
        a = parameters(1:m);
        b = parameters(m+1:m+n);
    case 'Diagonal'
        a = parameters(1:m*k);
        b = parameters(k*m+1:k*(m+n));
    case 'CP'
        a = parameters(1:m*k);
        lambda_cp=parameters(m*k+1);
        b= lambda_cp - a;
end

if stage==2
    intercept = R*(1-sum(a)-sum(b));
    [~,rescaledIntercept] = cov2corr(intercept);
    z = r2z(rescaledIntercept);
    startingVals = [z' parameters];
    
    LB = [-inf*ones(1,k*(k-1)/2) zeros(1,length(parameters))];
    UB = [inf*ones(1,k*(k-1)/2) ones(1,length(parameters))];
    A = [zeros(k,k*(k-1)/2) A];
    b = repmat(0.99998, k, 1);;
    parameters = fmincon(@rdcc_likelihood,startingVals,A,b,[],[],LB,UB,[],options,stdData,m,n,R,backCast,2,composite,isJoint,isInference,[],specification);
    z = parameters(1:k*(k-1)/2);
    R = z2r(z);
    parameters = parameters(k*(k-1)/2+1:length(parameters));
    parameters  = [corr_vech(R)' parameters];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariances and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
garchParameters = [];
for i=1:k
    garchParameters  = [garchParameters univariate{i}.parameters']; %#ok<AGROW>
end
if stage==3 
    parameters = [garchParameters corr_vech(R)' parameters];
elseif stage==2
    parameters = [garchParameters parameters];
end

isJoint = true; % parameters includes volatility parameters
isInference = true;
[ll,~,Rt] = rdcc_likelihood(parameters,data,m,n,R,backCast,stage,composite,isJoint,isInference,univariate, specification);
ll = -ll;
Ht = zeros(k,k,T);
for t=1:T
    h = sqrt(H(t,:));
    Ht(:,:,t) = Rt(:,:,t).*(h'*h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout<=3
    return
end

v = length(parameters);
A = zeros(v);
scores = zeros(T,v);
offset = 0;

for i=1:k
    u = univariate{i};
    count = 1 + u.p  + u.q;
    ind = offset+(1:count);
    A(ind,ind) = u.A;
    offset = offset + count;
    scores(:,ind) = u.scores;
end

% TODO : Better gradient function
l=0;
if stage==2
    % 1. dcc_likelihood
    
    count = k*(k-1)/2 + m  + n;
    H = hessian_2sided_nrows(@rdcc_likelihood,parameters',count,data,m,n,R,backCast,stage,composite,isJoint,isInference,univariate,specification);
    A(offset+(1:count),:) = H/T;
    [~,s]=gradient_2sided(@rdcc_likelihood,parameters',data,m,n,R,backCast,stage,composite,isJoint,isInference,univariate,specification);
    scores(:,offset+(1:count)) = s(:,offset+(1:count));
    B = cov(scores);
    Ainv = A\eye(v);
    VCV = Ainv*B*Ainv'/T;
elseif stage==3
    % 1. dcc_inference_objective
    count = k*(k-1)/2;
    tempParams = parameters(1:offset+count);
    DataAsym = zeros(k,k,T);
    [~,s]=gradient_2sided(@dcc_inference_objective, tempParams', data,DataAsym,m,l,n,univariate);
    scores(:,offset+(1:count))=s(:,offset+(1:count));
    H = hessian_2sided_nrows(@dcc_inference_objective, tempParams', count, data,DataAsym,m,l,n,univariate);
    A(offset+(1:count),1:(count+offset)) = H/T;
    offset = offset + count;
    % 2. dcc_likelihood
    count = m + l + n;
    H = hessian_2sided_nrows(@rdcc_likelihood,parameters',count,data,m,n,R,backCast,stage,composite,isJoint,isInference,univariate,specification);
    A(offset+(1:count),:) = H/T;
    [~,s]=gradient_2sided(@rdcc_likelihood,parameters',data,m,n,R,backCast,stage,composite,isJoint,isInference,univariate,specification);
    scores(:,offset+(1:count)) = s(:,offset+(1:count));
    B = covnw(scores);
    Ainv = A\eye(v);
    VCV = Ainv*B*Ainv'/T;
end


