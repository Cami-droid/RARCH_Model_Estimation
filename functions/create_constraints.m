function [A, b] = create_constraints(k, m, n)
    % k: number of assets
    % m: number of ARCH terms
    % n: number of GARCH terms
    
    % Total number of parameters per asset
    total_params_per_asset = m + n;
    
    % Total number of parameters for all assets
    total_params = k * total_params_per_asset;
    
    % Initialize A matrix and b vector
    A = zeros(k, total_params);
    b = ones(k, 1);  % We want the sum to be equal to 1 for each asset
    
    % Construct A matrix
    for i = 1:k
        % Indices for the parameters of the i-th asset
        start_idx = (i - 1) * total_params_per_asset + 1;
        end_idx_a = start_idx + m - 1;
        end_idx_b = start_idx + total_params_per_asset - 1;
        
        % Set the coefficients in the A matrix
        A(i, start_idx:end_idx_a) = 1;   % Coefficients for 'a' parameters
        A(i, end_idx_a+1:end_idx_b) = 1; % Coefficients for 'b' parameters
    end
end
