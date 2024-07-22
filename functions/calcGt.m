function Gt = calcGt(thetaS, thetaD, rotated_returns, Gt_prev, model, specification)
    %% thetaS: Static parameters. In this case, it's H_bar the unconditional covariance matrix (I don't define thetaS as vech(H_bar)) %%
    % thetaD: Dynamic parameters vector
    % rotated_returns: is a 1xd vector with returns previously rotated at time t-1
    % Gt_prev: G(t-1)
    % model: can be 'RBEKK', 'OGARCH', 'GOGARCH' or 'RDCC'
    % specification: 'Scalar', 'Diagonal', 'CP'

    d = size(rotated_returns, 2);
    
    % Extract dynamic parameters based on specification
    if strcmp(specification, 'Scalar')
        % thetaD = [alpha, beta] % 2 parameters
        alpha_vec = thetaD(1) * ones(d, 1); 
        beta_vec = thetaD(2) * ones(d, 1);
    elseif strcmp(specification, 'Diagonal')
        % thetaD = [alpha(1),...,alpha(d),beta(1),...,beta(d)] % 2d parameters
        alpha_vec = thetaD(1:d);
        beta_vec = thetaD(d+1:2*d);
    elseif strcmp(specification, 'CP')
        % thetaD = [alpha, lambda_cp] % d+1 parameters
        alpha_vec = thetaD(1) * ones(d, 1);
        beta_vec = (1 - thetaD(1)) * ones(d, 1);
    end
    
    A = diag(sqrt(alpha_vec));
    B = diag(sqrt(beta_vec));
    
    % Small regularization term to avoid singularity
    reg_term = 1e-6 * eye(d);

    % Compute Gt based on the model
    if strcmp(model, 'RBEKK')
        C = thetaS - A * thetaS * A' - B * thetaS * B';
        Ht = C + A * (rotated_returns' * rotated_returns) * A' + B * Gt_prev * B';
        Gt = Ht + reg_term; % Adding regularization term
    elseif strcmp(model, 'GOGARCH') || strcmp(model, 'OGARCH')
        C = eye(d) - A * A' - B * B';
        Gt = C + A * (rotated_returns' * rotated_returns) * A' + B * Gt_prev * B' + reg_term; % Adding regularization term
    elseif strcmp(model, 'RDCC')
        % Additional computations for RDCC if necessary  
        C = eye(d) - A * A' - B * B';
        Gt = C + A * (rotated_returns' * rotated_returns) * A' + B * Gt_prev * B' + reg_term; % Adding regularization term
    end
end
