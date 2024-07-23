function Gt = calcGt(model, specification, outputs, thetaD, Gt_prev)

    et = outputs.rotated_returns;
    thetaS = outputs.H_bar;
    d = outputs.d;

    % Extract dynamic parameters based on specification
    if strcmp(specification, 'Scalar')

        alpha_vec = thetaD(1) * ones(d, 1); 
        beta_vec = thetaD(2) * ones(d, 1);

    elseif strcmp(specification, 'Diagonal')

        alpha_vec = thetaD(1:d);
        beta_vec = thetaD(d+1:2*d);

    elseif strcmp(specification, 'CP')

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
        Ht = C + A * (et' * et) * A' + B * Gt_prev * B';

        Gt = Ht + reg_term; % Adding regularization term

    elseif strcmp(model, 'GOGARCH') || strcmp(model, 'OGARCH')

        C = eye(d) - A * A' - B * B';
        Gt = C + A * (et' * et) * A' + B * Gt_prev * B' + reg_term; % Adding regularization term

    elseif strcmp(model, 'RDCC')

        C = eye(d) - A * A' - B * B';
        Gt = C + A * (et' * et) * A' + B * Gt_prev * B' + reg_term; % Adding regularization term
    end

    % Check for NaN values and correct them
    if any(isnan(Gt(:)))
        warning('NaN values detected in Gt, adding additional regularization.');
        Gt = Gt + 1e-4 * eye(d);
    end

    % Ensure Gt is symmetric
    Gt = (Gt + Gt') / 2;
end
