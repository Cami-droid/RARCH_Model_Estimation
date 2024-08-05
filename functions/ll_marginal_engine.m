function L_marginal = ll_marginal_engine(model, outputs, thetaD)
    % ll_marginal_engine: Calculates the marginal log-likelihoods for each asset
    % Input arguments:
    %   model: String indicating the model type (e.g., 'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC')
    %   outputs: Struct with input data (r_t, etc.)
    %   thetaD: Dynamic parameter vector
    % Output arguments:
    %   L_marginal: Txd size matrix with marginal log-likelihoods

    % Initial variables
    rt = outputs.returns; % Centered return,  Txd matrix
    et = outputs.rotated_returns; % Rotated returns,  Txd matrix
    Et=outputs.std_returns; %standarized returns Txd matrix
    T = outputs.T; % Number of observations
    d = outputs.d; % Number of asset
    Dt = outputs.Dt; % Diagonal Matrices  with standard deviations , matrix dxdxT
    Lambda = outputs.Lambda; % eigenvalues for OGARCH y GOGARCH
    P = outputs.P; % eigenvectors OGARCH y GOGARCH

    % delta parameter for  OGARCH y GOGARCH
    delta = thetaD(end); 

    % Inicialize output matrices
    L_marginal = zeros(T, d);

    % calculate the marginal log-likelihood values
    switch model
        case 'RBEKK'
            for i = 1:d
                for t = 1:T
                    sigma2_t = outputs.Gt(i, i, t); % Conditional variance
                    L_marginal(t, i) = -0.5 * (log(2 * pi) + log(sigma2_t) + (et(t, i)^2) / sigma2_t);
                end
            end

        case {'OGARCH', 'GOGARCH'}
            U = @(delta, d) delta * eye(d);
            for i = 1:d
                for t = 1:T
                    sigma2_t = U(delta, d) * outputs.Gt(:,:,t) * U(delta, d)';
                    L_marginal(t, i) = -0.5 * (log(2 * pi) + log(sigma2_t(i, i)) + (et(t, i)^2) / sigma2_t(i, i));
                end
            end

        case 'RDCC'
            for i = 1:d
                for t = 1:T
                    sigma2_t = Dt(i, i, t)^2; % diagonal conditional variance
                    L_marginal(t, i) = -0.5 * (log(2 * pi) + log(sigma2_t) + (rt(t, i)^2) / sigma2_t);
                end
            end
    end
end 