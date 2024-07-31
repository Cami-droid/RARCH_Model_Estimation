function ll_Ms = ll_Ms_garch(returns)
    % data:  T x d de return matrix, where T is the number of observations and d is the number of series
    [T, d] = size(returns);
    ll_Ms = zeros(d, 1);

    for i = 1:d
        series = returns(:, i); % get i-esm series
        model = garch(1, 1); % define a  GARCH(1,1) model
        [~, ~, logL] = estimate(model, series, 'Display', 'off'); % Estimar el modelo y obtener la log-verosimilitud
        ll_Ms(i) = logL; % save the marginal log-likelihood
    end
end
