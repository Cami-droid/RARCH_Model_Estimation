function [parameters, ll ,Ht, VCV, scores, diagnostics] = rdcc(data, m, n, p, q, method, composite, specification, startingVals, options)
    % Estimation of the DCC(m,n) multivariate model with GARCH(p,q) conditional variances
    %
    % USAGE:
    %   [PARAMETERS] = dcc(DATA, M, N)
    %   [PARAMETERS, LL, HT, VCV, SCORES, DIAGNOSTICS] = rdcc(DATA, M, N, P, Q, METHOD, COMPOSITE, SPECIFICATION, STARTINGVALS, OPTIONS)
    %
    % INPUTS:
    %   DATA         - T by K matrix of zero mean residuals.
    %   M            - Order of symmetric innovations in the DCC model.
    %   N            - Order of lagged correlations in the DCC model.
    %   P            - Number of symmetric innovations in univariate volatility models.
    %   Q            - Number of conditional covariance lags in univariate volatility models.
    %   METHOD       - Estimation method ('3-stage' or '2-stage').
    %   COMPOSITE    - Composite likelihood type ('None', 'Diagonal', 'Full').
    %   SPECIFICATION- Model specification ('Scalar', 'Diagonal', 'CP').
    %   STARTINGVALS - Initial values.
    %   OPTIONS      - Optimization options (fmincon).
    %
    % OUTPUTS:
    %   PARAMETERS   - Estimated parameters.
    %   LL           - Log-likelihood at the optimum.
    %   HT           - Conditional covariance matrix [K K T].
    %   VCV          - Robust parameter covariance matrix.
    %   SCORES       - T by numParams matrix of individual scores.
    %   DIAGNOSTICS  - Structure with additional outputs.

    % Input validation
    if nargin < 3
        error('At least three inputs are required.');
    end
    if nargin < 4, p = 1; end
    if nargin < 5, q = 1; end
    if nargin < 6, method = '3-stage'; end
    if nargin < 7, composite = 'none'; end
    if nargin < 8, specification = 'Scalar'; end
    if nargin < 9, startingVals = []; end
    if nargin < 10, options = optimset('fmincon'); options.Display = 'iter'; options.Diagnostics = 'on'; options.Algorithm = 'interior-point'; end

    [T, K] = size(data);
    data2d = data;
    data = zeros(K, K, T);

    for t = 1:T
        data(:, :, t) = data2d(t, :)' * data2d(t, :);
    end

    % Parameter validation
    if ~isscalar(m) || floor(m) ~= m || m < 1
        error('M must be a positive integer.');
    end
    if ~isscalar(n) || floor(n) ~= n || n < 0
        error('N must be a non-negative integer.');
    end
    if ~isscalar(p), p = ones(K, 1); end
    if ~isscalar(q), q = ones(K, 1); end

    method = lower(method);
    if ~ismember(method, {'3-stage', '2-stage'})
        error('METHOD must be ''3-stage'' or ''2-stage''.');
    end
    stage = strcmp(method, '3-stage') + 1;

    composite = lower(composite);
    if ~ismember(composite, {'none', 'diagonal', 'full'})
        error('COMPOSITE should be ''none'', ''diagonal'' or ''full''.');
    end
    composite = find(strcmp(composite, {'none', 'diagonal', 'full'})) - 1;
    if stage == 2 && composite == 1
        warning('when METHOD is ''2-stage'', COMPOSITE should be ''none'' or ''full''.');
        composite = 2;
    end

    % Setting up initial values
    if isempty(startingVals)
        count = K + sum(p) + sum(q);
        count = count + m + n;
        if stage == 2
            count = count + K * (K - 1) / 2;
        end
        startingVals = rand(count, 1); % Replace this with more meaningful starting values if needed
    end

    % Univariate volatility models
    [H, univariate] = dcc_fit_variance(data2d, p, q, startingVals);
    stdData = data;

    for t = 1:T
        h = sqrt(H(t, :));
        hh = h' * h;
        stdData(:, :, t) = stdData(:, :, t) ./ hh;
    end

    % Estimation
    R = mean(stdData, 3);
    r = sqrt(diag(R));
    R = R ./ (r * r');

    % Setting up initial parameters for DCC
    if isempty(startingVals)
        a = [.01 .03 .05 .1];
        theta = [.99 .97 .95];
        [a, theta] = ndgrid(a, theta);
        parameters = unique([a(:) theta(:) - a(:)], 'rows');
        minLL = inf;

        for i = 1:size(parameters, 1)
            ll = dcc_spec_likelihood(parameters(i, :), stdData, [], m, l, n, R, [], [], [], stage, composite, false, false, [], univariate, specification);
            if ll < minLL
                startingVals = parameters(i, :);
                minLL = ll;
            end
        end
        a = startingVals(1);
        b = startingVals(2);
        startingVals = [a * ones(1, m) / m b * ones(1, n) / n];
    end

    LB = zeros(length(startingVals), 1);
    UB = ones(length(startingVals), 1);
    A = ones(1, length(startingVals));
    b = .99998;

    if (startingVals * A' - b) >= 0
        error('STARTINGVALS for DCC parameters are not compatible with a positive definite intercept.');
    end

    % Optimization
    parameters = fmincon(@(params) dcc_spec_likelihood(params, stdData, [], m, l, n, R, [], [], [], stage, composite, false, false, [], univariate, specification), startingVals, A, b, [], [], LB, UB, [], options);

    % Covariances and parameters
    garchParameters = [];
    for i = 1:K
        garchParameters = [garchParameters univariate{i}.parameters']; %#ok<AGROW>
    end

    if stage == 2
        intercept = R * (1 - sum(parameters(1:m)) - sum(parameters(m+1:m+n)));
        [~, rescaledIntercept] = cov2corr(intercept);
        z = r2z(rescaledIntercept);
        startingVals = [z' parameters];

        LB = [-inf * ones(1, K * (K - 1) / 2) zeros(1, length(parameters))];
        UB = [inf * ones(1, K * (K - 1) / 2) ones(1, length(parameters))];
        A = [zeros(1, K * (K - 1) / 2) A];
        parameters = fmincon(@(params) dcc_spec_likelihood(params, stdData, [], m, l, n, R, [], [], [], stage, composite, false, false, [], univariate, specification), startingVals, A, b, [], [], LB, UB, [], options);
        z = parameters(1:K * (K - 1) / 2);
        R = z2r(z);
        parameters = parameters(K * (K - 1) / 2 + 1:end);
        parameters = [corr_vech(R)' parameters];
    end

    parameters = [garchParameters parameters];

    [ll, ~, Rt] = dcc_spec_likelihood(parameters, data, [], m, l, n, R, [], [], [], stage, composite, true, true, [], univariate, specification);
    ll = -ll;
    Ht = zeros(K, K, T);

    for t = 1:T
        h = sqrt(H(t, :));
        Ht(:, :, t) = Rt(:, :, t) .* (h' * h);
    end

    % Inference
    if nargout > 3
        v = length(parameters);
        A = zeros(v);
        scores = zeros(T, v);
        offset = 0;

        for i = 1:K
            u = univariate{i};
            count = 1 + u.p + u.q;
            ind = offset + (1:count);
            A(ind, ind) = u.A;
            offset = offset + count;
            scores(:, ind) = u.scores;
        end

        count = K * (K - 1) / 2 + m + n;
        H = hessian_2sided_nrows(@(params) dcc_spec_likelihood(params, data, [], m, l, n, R, [], [], [], stage, composite, true, true, [], univariate, specification), parameters', count, data, m, n, R, stage, composite, true, true);
        A(offset + (1:count), :) = H / T;
        [~, s] = gradient_2sided(@(params) dcc_spec_likelihood(params, data, [], m, l, n, R, [], [], [], stage, composite, true, true, [], univariate, specification), parameters', data, m, n, R, stage, composite, true, true);
        scores(:, offset + (1:count)) = s(:, offset + (1:count));
        B = cov(scores);
        Ainv = A \ eye(v);
        VCV = Ainv * B * Ainv' / T;
    end

    diagnostics = [];
end

% Required functions as  dcc_fit_variance, dcc_spec_likelihood, hessian_2sided_nrows, gradient_2sided,
% corr_vech, z2r, r2z, and others have to be defined separately or imported as needed.
