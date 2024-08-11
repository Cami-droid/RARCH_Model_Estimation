function [ll, lls, Rt] = dcc_spec_likelihood(parameters, data, dataAsym, m, l, n, R, N, backCast, backCastAsym, stage, composite, isJoint, isInference, gScale, univariate, specification)
    % Likelihood for estimation of scalar DCC(m,n) and ADCC(m,l,n) multivariate volatility models
    % with TARCH(p,o,q) or GJRGARCH(p,o,q) conditional variances
    %
    % USAGE:
    %  [LL,LLS,RT] = dcc_likelihood(PARAMETERS,DATA,DATAASYM,M,L,N,R,N,BACKCAST,BACKCASTASYM,STAGE,COMPOSITE,ISJOINT,ISINFERENCE,GSCALE,UNIVARIATE)
    %
    % INPUTS:
    %   PARAMETERS   - Vector of DCC parameters, and possibly volatility and intercept parameters,
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
    %   BACKCAST     - K by K  matrix to use for back casting symmetric terms
    %   BACKCASTASYM - K by K  matrix to use for back casting asymmetric terms
    %   STAGE        - Integer, either 2 or 3 indicating whether 2-stage or 3-stage estimator is being used
    %   COMPOSITE    - Integer, one of 0 (None, use QMLE), 1 (Use diagonal composite) or 2 (full composite)
    %   ISJOINT      - Boolean indicating whether PARAMETERS includes volatility parameters
    %   ISINFERENCE  - Boolean indicating whether likelihood is used for making inference, in which case
    %                    no transformations are made to parameters.
    %   GSCALE       - K by 1 vector used in 2-stage to scale the intercept.  See DCC.
    %   UNIVARIATE   - Cell array of structures containing information needed to compute volatilities.
    %   SPECIFICATION- String, one of 'Scalar', 'Diagonal' or 'CP'
    %
    % OUTPUTS:
    %   LL           - The log likelihood evaluated at the PARAMETERS
    %   LLS          - A T by 1 vector of log-likelihoods
    %   RT           - A [K K T] dimension matrix of conditional correlations
    %
    % COMMENTS:
    %
    % See also DCC

    [k,~,T] = size(data);
    offset = 0;

    % Parse Parameters
    if stage == 1 || isJoint
        count = 0;
        for i = 1:k
            u = univariate{i};
            count = count + u.p + u.o + u.q + 1;
        end
        garchParameters = parameters(1:count);
        offset = offset + count;
        computeVol = true;
    else
        computeVol = false;
    end

    if stage <= 2 || isJoint
        count = k * (k - 1) / 2;
        R = parameters(offset + (1:count));
        offset = offset + count;
    end

    if stage==3 && isJoint && l>0 % N only if 3 stage, joint and asymmetric
    count = k*(k+1)/2;
    N = parameters(offset + (1:count));
    % Transform N to be a matrix
    N = ivech(N);
    offset = offset + count;
    end

    switch specification
        case 'Scalar'
            a = parameters(offset + (1:m));
            b = parameters(offset + (m+1:m+n));
        case 'Diagonal'
            a = parameters(offset + (1:m*k));
            b = parameters(offset + (m*k+1:m*k+n*k));
            a = reshape(a, [k, m]);
            b = reshape(b, [k, n]);
        case 'CP'
            a = parameters(offset + (1:m*k));
            lambda_cp = parameters(offset + m*k + 1);
            b=(lambda_cp*ones(k,1)-a.^4).^(1/4);
            
    end

    % Compute volatilities
    H = ones(T, k);
    if computeVol
        H = dcc_reconstruct_variance(garchParameters, univariate);
        stdData = zeros(k, k, T);
        stdDataAsym = zeros(k, k, T);
        for t = 1:T
            h = sqrt(H(t, :));
            stdData(:, :, t) = data(:, :, t) ./ (h' * h);
            stdDataAsym(:, :, t) = dataAsym(:, :, t) ./ (h' * h);
        end
        logdetH = sum(log(H), 2);
    else
        stdData = data;
        stdDataAsym = dataAsym;
        logdetH = zeros(T, 1);
    end

    % Transform R if needed
    if stage <= 2 || isJoint
        if isInference
            R = corr_ivech(R);
        else
            R = z2r(R);
        end
    end

    % Compute intercept
    if stage == 3
        intercept = R * (1 - sum(a(:)) - sum(b(:)));
        if l > 0
            intercept = intercept - N * sum(g);
        end
    else
        if l > 0 && exist('g', 'var')
            scale = (1 - sum(a(:)) - sum(b(:))) - gScale * sum(g);
        else
            scale = 1 - sum(a(:)) - sum(b(:));
        end
        scale = sqrt(scale);
        intercept = R .* (scale * scale');
    end


    % Check eigenvalues?

    % Indices or constant, as needed
    if composite == 0
        likconst = k * log(2 * pi);
    elseif composite == 1
        indices = [(1:k-1)' (2:k)'];
    elseif composite == 2
        [i, j] = meshgrid(1:k);
        indices = [i(~triu(true(k))) j(~triu(true(k)))];
    end

    I = eye(k);
    Qt = zeros(k, k, T);
    Rt = zeros(k, k, T);
    lls = zeros(T, 1);

    for t = 1:T
        Qt(:, :, t) = intercept;
        for i = 1:m
            if (t-i) > 0
                if strcmp(specification, 'Scalar')
                    Qt(:, :, t) = Qt(:, :, t) + a(i) * stdData(:, :, t-i);
                else
                    Qt(:, :, t) = Qt(:, :, t) + diag(a(:, i)) * stdData(:, :, t-i) * diag(a(:, i));
                end
            else
                if strcmp(specification, 'Scalar')
                    Qt(:, :, t) = Qt(:, :, t) + a(i) * backCast;
                else
                    Qt(:, :, t) = Qt(:, :, t) + diag(a(:, i)) * backCast * diag(a(:, i));
                end
            end
        end

        if strcmp(specification, 'CP')
            for i = 1:n
                if (t-i) > 0
                    Qt(:, :, t) = Qt(:, :, t) + lambda_cp * Qt(:, :, t-i);
                else
                    Qt(:, :, t) = Qt(:, :, t) + lambda_cp * backCast;
                end
            end
        else
            for i = 1:n
                if (t-i) > 0
                    Qt(:, :, t) = Qt(:, :, t) + diag(b(:, i)) * Qt(:, :, t-i) * diag(b(:, i));
                else
                    Qt(:, :, t) = Qt(:, :, t) + diag(b(:, i)) * backCast * diag(b(:, i));
                end
            end
        end

        q = sqrt(diag(Qt(:, :, t)));
        Rt(:, :, t) = Qt(:, :, t) ./ (q * q');
        if composite == 0
            lls(t) = 0.5 * (likconst + logdetH(t) + log(det(Rt(:, :, t))) + sum(diag((Rt(:, :, t) \ I) * stdData(:, :, t))));
        elseif composite
            S = (sqrt(H(t, :))' * sqrt(H(t, :))) .* Rt(:, :, t);
            lls(t) = composite_likelihood(S, data(:, :, t), indices);
        end
    end

    ll = sum(lls);

    if isnan(ll) || ~isreal(ll) || isinf(ll)
        ll = 1e7;
    end
end

% Subfunctions required as dcc_fit_variance, dcc_likelihood, hessian_2sided_nrows, gradient_2sided,
% corr_vech, z2r, r2z, y otras se deben definir por separado o importar según sea necesario.
