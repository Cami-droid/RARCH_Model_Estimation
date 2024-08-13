function [parameters, ll ,Ht, VCV, scores, diagnostics]=rdcc(data,dataAsym,m,l,n,p,o,q,gjrType,method,composite,startingVals,options,specification)
% Estimation of scalar DCC(m,n) and ADCC(m,l,n) multivarate volatility model with with TARCH(p,o,q)
% or GJRGARCH(p,o,q) conditional variances
%
% USAGE:
%  [PARAMETERS] = dcc(DATA,[],M,L,N)
%  [PARAMETERS,LL,HT,VCV,SCORES,DIAGNOSTICS] =
%  dcc(DATA,DATAASYM,M,L,N,P,O,Q,
%                                                    GJRTYPE,METHOD,COMPOSITE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   DATAASYM     - [OPTIONAL] K by K by T array of asymmetric covariance estimators only needed if
%                    DATA is 3-dimensional and O>0 or L>0
%   M            - Order of symmetric innovations in DCC model
%   L            - Order of asymmetric innovations in ADCC model
%   N            - Order of lagged correlation in DCC model
%   P            - [OPTIONAL] Positive, scalar integer representing the number of symmetric innovations in the
%                    univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 1.
%   O            - [OPTIONAL] Non-negative, scalar integer representing the number of asymmetric innovations in the
%                    univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 0.
%   Q            - [OPTIONAL] Non-negative, scalar integer representing the number of conditional covariance lags in
%                    the univariate volatility models.  Can also be a K by 1 vector containing the lag length
%                    for each series. Default is 1.
%   GJRTYPE      - [OPTIONAL] Either 1 (TARCH/AVGARCH) or 2 (GJR-GARCH/GARCH/ARCH). Can also be a K by 1 vector
%                    containing the model type for each for each series. Default is 2.
%   METHOD       - [OPTIONAL] String, one of '3-stage' (Default) or '2-stage'.  Determines whether
%                    the model is estimated using the 3-stage estimator, or if the correlation intercepts
%                    are jointly estimated along with the dynamic parameters.
%   COMPOSITE    - [OPTIONAL] String value, either 'None' (Default), 'Diagonal' or 'Full'.  None
%                    uses standard QMLE.  'Diagonal' and 'Full' both uses composite likelihood where
%                    'Diagonal' uses all pairs of the form i,i+1 while 'Full' uses all pairs.
%   STARTINGVALS - [OPTIONAL] Vector of starting values to use.  See parameters and COMMENTS.
%   OPTIONS      - [OPTIONAL] Options to use in the model optimization (fmincon)
%   SPECIFICATION -[OPTIONAL] String value, either 'Scalar'(Default),'Diagonal','CP' (Common Persistence). When 'CP' is selected, l is set to zero 

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
%   where v(t,:) = e(t,;).*(e(t,:)<0) and s = sqrt((1-sum(a)-sum(b)-gScale*sum(g))) and scale = s*s'
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
    case 5
        p = [];
        o = [];
        q = [];
        gjrType = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
        specification=[];
    case 6
        o = [];
        q = [];
        gjrType = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
        specification=[];
    case 7
        q = [];
        gjrType = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
        specification=[];
    case 8
        gjrType = [];
        method = [];
        composite = [];
        startingVals = [];
        options = [];
        specification=[];
    case 9
        method = [];
        composite = [];
        startingVals = [];
        options = [];
        specification=[];
    case 10
        composite = [];
        startingVals = [];
        options = [];
        specification=[];
    case 11
        startingVals = [];
        options = [];
        specification=[];
    case 12
        options = [];
        specification=[];
    case 13
        specification=[];
    case 14
    % Nothing
    otherwise
        error('5 to 14 inputs required.')
end

if ismatrix(data)
    [T,k] = size(data);
    data2d = data;
    eta = data.*(data<0);
    data = zeros(k,k,T);
    dataAsym = zeros(k,k,T);
    for t=1:T
        data(:,:,t) = data2d(t,:)'*data2d(t,:);
        dataAsym(:,:,t) = eta(t,:)'*eta(t,:);
    end
elseif ndims(data)==3
    [k,~,T]=size(data);
    data2d = zeros(k,T);
    for i=1:k
        data2d(:,i) = squeeze((1-2*(dataAsym(i,i,:)==0)) .* sqrt(data(i,i,:)));
    end
else
    error('DATA must be either a K by T matrix of a K by K by T 3-dimensional array.')
end

if ~isempty(dataAsym)
    if ndims(dataAsym)~=3
        error('DATAASYM must be a K by K by T array.')
    end
    [k1,k2,T2] = size(dataAsym);
    if k1~=k || k2~=k || T2~=T
        error('DATAASYM must be a K by K by T array.')
    end
end

if ~isscalar(m) || floor(m)~=m || m<1
    error('M must be a positive integer.')
end
if isempty(l)
    l = 0;
end

% 'CP' is not developed for l>0 in this code. if its value is different than zero, automatically is set to be zero.
if strcmp(specification,'CP')
l=0; 
g=0;
end
if ~isscalar(l) || floor(l)~=l || l<0
    error('L must be a non-negative integer.')
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
if isempty(o)
    o = 0;
end
if isscalar(o)
    o = o * ones(k,1);
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
if any(floor(o)~=o) || any(o<0) || numel(o)~=k
    error('All elements in O must be non-negative, and O must be either scalar or K by 1.')
end
if any(floor(q)~=q) || any(q<0) || numel(q)~=k
    error('All elements in Q must be non-negative, and Q must be either scalar or K by 1.')
end
if isempty(gjrType)
    gjrType = 2;
end
if isscalar(gjrType)
    gjrType = ones(k,1)*gjrType;
end
if any(~ismember(gjrType,[1 2])) || numel(gjrType)~=k
    error('GJRTYPE must be in {1,2} and must be either scalar of K by 1.')
end

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
switch specification
    case 'Scalar'
    count_thetaD=m+l+n;
    case 'Diagonal'
    count_thetaD=k*(m+l+n);
    case 'CP'
    count_thetaD=k*m+1;
    otherwise
    count_thetaD=m+l+n;
end

if ~isempty(startingVals)
    count = k + sum(p)+ sum(o) + sum(q);
    count = count + m + l + n;
    if stage==2 || l==0
        count = count + k*(k-1)/2;
    elseif  stage==3
        count = count + k*(k+1)/2 + k*(k-1)/2;
    end
    if length(startingVals)~=count
        error('STARTINGVALS does not contain the correct number of parameters.')
    end
    count = k+sum(p)+sum(o)+sum(q);
    tarchStartingVals = startingVals(1:count);
    offset = count;
    % Count for intercepts
    if stage==2 || l==0
        count = k*(k-1)/2;
    elseif  stage==3
        count = k*(k+1)/2 + k*(k-1)/2;
    end
    offset = offset + count;
    count = count_thetaD;
    dccStartingVals = startingVals(offset + (1:count_thetaD));
else
    tarchStartingVals = [];
    dccStartingVals = [];
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

if isempty(specification);
specification='Scalar'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Univariate volatility models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[H,univariate] = dcc_fit_variance(data2d,p,o,q,gjrType,tarchStartingVals);
stdData = data;
stdDataAsym = dataAsym;
for t=1:T
    h = sqrt(H(t,:));
    hh = h'*h;
    stdData(:,:,t) = stdData(:,:,t)./hh;
    stdDataAsym(:,:,t) = stdDataAsym(:,:,t)./hh;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdData2d=data2d;
for t=1:T
    h=sqrt(H(t,:));
    stdData2d(t,:)=data2d(t,:)./h;  
end

PI_bar = mean(stdData,3);    % same value as stdData2d'*stdData2d/T
[P_C Lambda_C]=eig(PI_bar);
et2d=P_C*Lambda_C^(-0.5)*P_C'*stdData2d';
et2d=et2d';

etet=zeros(k,k,T);
for t=1:T
        etet(:,:,t) = et2d(t,:)'*et2d(t,:);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Back casts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = .06*.94.^(0:sqrt(T));
w = w'/sum(w);
backCast = zeros(k);
backCastAsym = zeros(k);
for i=1:length(w)
    backCast = backCast + w(i)*etet(:,:,i);
    backCastAsym = backCastAsym + w(i)*stdDataAsym(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = PI_bar;
r = sqrt(diag(R));
R = R ./ (r*r');
N = mean(stdDataAsym,3);
scale = max(eig(R^(-0.5)*N*R^(-0.5)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(dccStartingVals)
    startingVals = dccStartingVals;
else
    a = [.01 .03 .05 .1];
    theta = [.99 .97 .95];
    if l>0
        g = [.01 .03 .05];
        [a,g,theta] = ndgrid(a,g,theta);
        parameters = unique([a(:) g(:) theta(:)-a(:)-scale*g(:)],'rows');
        isAsym = 1;
    else
        [a,theta] = ndgrid(a,theta);
        parameters = unique([a(:) theta(:)-a(:)],'rows');
        switch specification
            case 'Scalar'
                % nothing to do here, parameters doesn't change
            case 'Diagonal'
                parameters1= repmat(parameters(:,1),1,k);
                parameters2= repmat(parameters(:,2),1,k);
                parameters=horzcat(parameters1,parameters2);
            clear parameters1;
            clear parameters2;
            case 'CP'
                parameters1= repmat(parameters(:,1),1,k);
                parameters2=parameters(:,2);
                parameters=horzcat(parameters1,parameters2);
                clear parameters1;
                clear parameters2;
        end 
        isAsym = 0;
    end
    minLL= inf;
    isJoint = false;
    isInference = false;
    for i=1:size(parameters,1)
        ll = dcc_spec_likelihood(parameters(i,:),etet,stdDataAsym,1,isAsym,1,R,N,backCast,backCastAsym,3,composite,isJoint,isInference,[],options,specification);
        if ll<minLL
            startingVals = parameters(i,:);            
            minLL = ll;
        end
    end
    switch specification
        case 'Scalar'
            a = startingVals(1);
            if l>0
                g = startingVals(2);
                b = startingVals(3);
                startingVals = [a*ones(1,m)/m g*ones(1,l)/l b*ones(1,n)/n];
            else
                b = startingVals(2);
                startingVals = [a*ones(1,m)/m b*ones(1,n)/n];
            end
        case 'Diagonal'
            a = startingVals(1:m*k);
                if l>0
                    g = startingVals(m*k+1:(m+l)*k);
                    b = startingVals((m+l)*k+1:(m+l+n)*k);
                    startingVals = [a g b];
                else
                   b = startingVals((m+l)*k+1:(m+l+n)*k);
                    startingVals = [a b];
                end
        case 'CP'
            a= startingVals(1:m*k);
            lambda_cp= startingVals( m*k + 1);
            startingVals = [a lambda_cp];
        

    end


    end
    switch specification
        case 'Scalar'
            A(m+1:m+l) = scale;  % Applied to g
        case 'Diagonal'
            A((m*k+1):(m*k+l*k)) = scale;  % Applied to g
        case 'CP'
        % A debe reflejar la relación lambda_cp - a
            A((m*k+1):(m*k+l*k)) = scale;
    end
    b = 0.99998;

    LB = zeros(length(startingVals),1);
    UB = ones(length(startingVals),1);
    A = ones(1,length(startingVals));
    A(m+1:m+l) = scale;
    b = .99998;

    if (startingVals*A'-b) >= 0
        coef=(startingVals*A'-b)./b;
        startingVals=startingVals./coef; %to no obstaculize;
        % error('STARTINGVALS for DCC parameters are not compatible with a positive definite intercept.')
        clear coef;
    end

    isJoint = false;
    isInference = false;
    parameters = fmincon(@dcc_spec_likelihood,startingVals,A,b,[],[],LB,UB,[],options,etet,stdDataAsym,m,l,n,R,N,backCast,backCastAsym,3,composite,isJoint,isInference,[],options,specification);

    gScale = diag(N);

    switch specification
    case 'Scalar'
        a = parameters(1:m);
        g = parameters(m+1:m+l);
        b = parameters(m+l+1:m+l+n);
    case 'Diagonal'
        a = parameters(1:m*k);
        g = parameters(m*k+1:(m+l)*k);    
        b = parameters((m+l)*k+1:(m+l+n)*k);
    case 'CP'
        a = parameters(1:m*k);
        lambda_cp = parameters(m*k+1);
        b= lambda_cp-a;
    otherwise
        a = parameters(1:m);
        g = parameters(m+1:m+l);
        b = parameters(m+l+1:m+l+n);
end

if stage==2
    intercept = R*(1-sum(a)-sum(b)) - N*sum(g);
    [~,rescaledIntercept] = cov2corr(intercept);
    z = r2z(rescaledIntercept);
    startingVals = [z' parameters];
    
    LB = [-inf*ones(1,k*(k-1)/2) zeros(1,length(parameters))];
    UB = [inf*ones(1,k*(k-1)/2) ones(1,length(parameters))];
    A = [zeros(1,k*(k-1)/2) A];
    b = .99998;
    parameters = fmincon(@dcc_spec_likelihood,startingVals,A,b,[],[],LB,UB,[],options,etet,stdDataAsym,m,l,n,R,N,backCast,backCastAsym,2,composite,isJoint,isInference,gScale,[],specification);
    
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
    if l>0
        parameters = [garchParameters corr_vech(R)' vech(N)' parameters];
    else
        parameters = [garchParameters corr_vech(R)' parameters];
    end
elseif stage==2
    parameters = [garchParameters parameters];
end

isJoint = true;
isInference = true;
[ll,~,Rt] = dcc_spec_likelihood(parameters,data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate,specification);
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
    count = 1 + u.p + u.o + u.q;
    ind = offset+(1:count);
    A(ind,ind) = u.A;
    offset = offset + count;
    scores(:,ind) = u.scores;
end

% TODO : Better gradient function
if stage==2
    % 1. dcc_likelihood
    count = k*(k-1)/2 + m + l + n;
    H = hessian_2sided_nrows(@dcc_spec_likelihood,parameters',count,data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate,specification);
    A(offset+(1:count),:) = H/T;
    [~,s]=gradient_2sided(@dcc_spec_likelihood,parameters',data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate,specification);
    scores(:,offset+(1:count)) = s(:,offset+(1:count));
    B = cov(scores);
    Ainv = A\eye(v);
    VCV = Ainv*B*Ainv'/T;
elseif stage==3
    % 1. dcc_inference_objective
    count = k*(k-1)/2;
    if l>0
        count = count + k*(k+1)/2;
    end
    tempParams = parameters(1:offset+count);
    [~,s]=gradient_2sided(@dcc_inference_objective, tempParams', data,dataAsym,m,l,n,univariate);
    scores(:,offset+(1:count))=s(:,offset+(1:count));
    H = hessian_2sided_nrows(@dcc_inference_objective, tempParams', count, data,dataAsym,m,l,n,univariate);
    A(offset+(1:count),1:(count+offset)) = H/T;
    offset = offset + count;
    % 2. dcc_likelihood
    count = m + l + n;
    H = hessian_2sided_nrows(@dcc_likelihood,parameters',count,data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate);
    A(offset+(1:count),:) = H/T;
    [~,s]=gradient_2sided(@dcc_likelihood,parameters',data,dataAsym,m,l,n,[],[],backCast,backCastAsym,stage,composite,isJoint,isInference,gScale,univariate);
    scores(:,offset+(1:count)) = s(:,offset+(1:count));
    B = covnw(scores);
    Ainv = A\eye(v);
    VCV = Ainv*B*Ainv'/T;
end

diagnostics = [];
diagnostics.gScale = gScale;