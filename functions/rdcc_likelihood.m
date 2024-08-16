function [ll,lls,Rt] = rdcc_likelihood(parameters,data,m,n,R,backCast,stage,composite,isJoint,isInference,univariate,specification)
% Likelihood for estimation of scalar DCC(m,n) and ADCC(m,l,n) multivarate volatility models 
% with with TARCH(p,o,q) or GJRGARCH(p,o,q) conditional variances
%
% USAGE:
%  [LL,LLS,RT] = dcc_likelihood(PARAMETERS,DATA,M,L,N,R,N,BACKCAST,STAGE,COMPOSITE,ISJOINT,ISINFERENCE,UNIVARIATE)
%
% INPUTS:
%   PARAMETERS   - Vector of ADCC parameters, and possibly volatility and intercept parameters, 
%                    depending on other inputs
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   DATAASYM     - [OPTIONAL] K by K by T array of asymmetric covariance estimators only needed if
%                    DATA is 3-dimensional and O>0 or L>0
%   M            - Order of symmetric innovations in DCC model
%   L            - Order of asymmetric innovations in ADCC model
%   N            - Order of lagged correlation in DCC model
%   R            - K by K correlation matrix of standardized data
%   N            - K by K matrix mean of asymmetric standardized data outer products
%   BACKCAST     - K by K  matrix to use for back casting symetric terms
%   STAGE        - Integer, either 2 or 3 indicating whether 2-stage ro 3-stage estimator is being used
%   COMPOSITE    - Integer, one of 0 (None, use QMLE), 1 (Use diagonal composite) or 2 (full composite)
%   ISJOINT      - Boolean indicating whether PARAMETERS includes volatility parameters
%   ISINFERENCE  - Boolean indicating whether likelihood is used for making inference, in which case
%                    no transformations are made to parameters.
%   UNIVARIATE   - Cell array of structures containing information needed to compute volatilities.
%
% OUTPUTS:
%   LL           - The log likelihood evaluated at the PARAMETERS
%   LLS          - A T by 1 vector of log-likelihoods
%   RT           - A [K K T] dimension matrix of conditional correlations
%
% COMMENTS:
%
% See also DCC

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/13/2012


% Modes:
% 1-stage: Univariate, R (corr vech), DCC (N-scale ignored)
% 2-stage (Estimation): R (transformed corr vech), DCC  (N-scale treated as constants)
% 2-stage (Inference): Univariate, R (corr vech), DCC (N-scale ignored)
% 3-stage (Estimation): DCC
% 3-stage (Inferecne): Univariate, R (corr vech), N (vech), DCC

[k,~,T] = size(data);
offset = 0;
% Parse Parameters
if stage==1 || isJoint  % include marginal parameters
    count = 0;
    for i=1:k
        u = univariate{i};
        count = count + u.p+u.q+1;
    end
    garchParameters = parameters(1:count);
    offset = offset + count;
    computeVol = true;
else
    computeVol = false;
end
if stage<=2 || isJoint
    count = k*(k-1)/2;
    R = parameters(offset + (1:count));
    offset = offset + count;
end

switch specification

    case 'Scalar'
        a = parameters((1:m));
        b = parameters(m+1:m+n);

          
    case 'Diagonal'
        a = parameters(1:m*k);
        b = parameters(k*m+1:k*(m+n));

        a=reshape(a,k,[]);
        b=reshape(b,k,[]);

    case 'CP'
        
        a = parameters(1:m*k);
        lambda_cp=parameters(m*k+1);
        b=sqrt(lambda_cp-a.^2);

        a=reshape(a,k,[]);
        b=reshape(b,k,[]);

end

if isempty(b)
    b=0;
end
% Compute volatilities
H = ones(T,k);
if computeVol
    H = dcc_reconstruct_variance(garchParameters,univariate);
    stdData = zeros(k,k,T);
    for t=1:T
        h = sqrt(H(t,:));
        stdData(:,:,t) = data(:,:,t)./(h'*h);
    end
    logdetH = sum(log(H),2);
else
    stdData = data;
    logdetH = zeros(T,1);
end




% Transfor R & N, if needed
if stage<=2 || isJoint
    % Transform R
    if isInference
        R = corr_ivech(R);
    else
        R = z2r(R);
    end
end

% Compute intercept
if stage==3
    intercept = R*(1-sum(a)-sum(b));
else
    scale = (1-sum(a)-sum(b));
    scale = sqrt(scale);
    intercept = R.*(scale*scale');
end
% Check eigenvalues?

% Indices or constant, as needed
if composite == 0
    likconst = k*log(2*pi);
elseif composite == 1
    indices = [(1:k-1)' (2:k)'];
elseif composite == 2
    [i,j] = meshgrid(1:k);
    indices = [i(~triu(true(k))) j(~triu(true(k)))];
end

I = eye(k);
Qt = zeros(k,k,T);
Rt = zeros(k,k,T);
lls = zeros(T,1);
for t=1:T
    Qt(:,:,t) = intercept;
    for i = 1:m
        if (t-i)>0
            Qt(:,:,t) = Qt(:,:,t) +diag(a(:,i))*stdData(:,:,t-i)*diag(a(:,i))';
        else
            Qt(:,:,t) = Qt(:,:,t) + diag(a(:,i))*backCast*diag(a(:,i))';
        end
    end
 
    for i = 1:n
        if (t-i)>0
            Qt(:,:,t) = Qt(:,:,t) + diag(b(:,i))*Qt(:,:,t-i)*diag(b(:,i))';
        else
            Qt(:,:,t) = Qt(:,:,t) + diag(b(:,i))*backCast*diag(b(:,i))';
        end
    end
    q = sqrt(diag(Qt(:,:,t)));
    Rt(:,:,t) = Qt(:,:,t)./ (q*q');
    if composite == 0
        lls(t) = 0.5*(likconst + logdetH(t) + log(det(Rt(:,:,t))) + sum(diag((Rt(:,:,t)\I)*stdData(:,:,t))));
    elseif composite
        S = (sqrt(H(t,:))'*sqrt(H(t,:))) .* Rt(:,:,t);
        lls(t) = composite_likelihood(S,data(:,:,t),indices);
    end
end
ll = sum(lls);

if isnan(ll) || ~isreal(ll) || isinf(ll)
    ll = 1e7;
end