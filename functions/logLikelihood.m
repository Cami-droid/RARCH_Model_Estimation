function ll = logLikelihoodFunc(thetaS, thetaD, rotated_returns, model, specification)
    T = size(rotated_returns, 1);
    d = size(rotated_returns, 2);
    
    % Preallocate Gt for efficiency
    Gt = zeros(d, d, T);
    Gt(:, :, 1) = eye(d); % Initial value of Gt for t=0
    
    ll = 0;
    % Gt = calcGt(thetaS,thetaD, rotated_returns, Gt_prev, model, specification)
    for t = 2:T+1 % Adjusting loop to account for t starting from 0
%         Gt(:, :, t) = calcGt(thetaS, thetaD, rotated_returns(t-1, :)', Gt(:, :, t-1), model, specification); % t-1 for returns to match theory
        if strcmp(model, 'GOGARCH')
            U_delta_et = rotated_returns(t-1, :)'; % Reemplazar con la operación específica U(δ)et para GOGARCH
            ll = ll - 0.5 * (d * log(2*pi) + log(det(Gt(:, :, t))) + U_delta_et * inv(Gt(:, :, t)) * U_delta_et');
        else
            ll = ll - 0.5 * (d * log(2*pi) + log(det(Gt(:, :, t))) + rotated_returns(t-1, :) * inv(Gt(:, :, t)) * rotated_returns(t-1, :)');
    end
end
