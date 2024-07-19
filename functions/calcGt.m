function Gt = calcGt(thetaS,thetaD, rotated_returns, Gt_prev, model, specification)
%   thetaS Static parameters. In this case is H_bar the unconditional covariance matrix
%%  thetaD dynamic parameters vector 

    % Scalar                thetaD=(alpha_vec beta_vec)
    % Diagonal              thetaD=(alpha_vec beta_vec)
    % Common Persistence    thetaD=(alpha_vec lambda_cp)
%%
% rotated_returns is a Txd matrix with returns previously rotated
% Gt_prev is G(t-1)
% model can be 'RBEKK', 'OGARCH', 'GOGARCH' AND 'RDCC'
% specification is 'Scalar', 'Diagonal', 'CP


    d = size(rotated_returns,2);
    
    % Extract dynamic parameters based on specification
    if strcmp(specification, 'Scalar')
        % thetaD=(alpha beta) % 2 parameters
        alpha_vec = thetaD(1) * ones(d, 1); 
        beta_vec = thetaD(2) * ones(d, 1);
    elseif strcmp(specification, 'Diagonal')
        % thetaD=(alpha(1)...(alpha(d),beta(1)...beta(d)) %2d parameters
        alpha_vec = thetaD(1:d);
        beta_vec = thetaD(d+1:2*d);
    elseif strcmp(specification, 'CP')
        % thetaD=(alpha(1)...(alpha(d),lambda_cp)   % d+1 parameters
        alpha_vec = thetaD(1) * ones(d, 1);
        beta_vec = (1 - thetaD(1)) * ones(d, 1);
    end
    
    A = diag(sqrt(alpha_vec));
    B = diag(sqrt(beta_vec));
    
    % Compute Gt based on the model
    if strcmp(model, 'RBEKK') 
        C = thetaS - A *thetaS* A' - B *thetaS* B';
        Ht = C + A * (rotated_returns' * rotated_returns) * A' + B * Gt_prev * B';
        Gt=Ht % this is not like that in theory, just to order the program that Ht is the output in the RBEKK case
        
        %FALTA COMPLETAR
    elseif strcmp(model, 'GOGARCH')|| strcmp(model, 'OGARCH')
        C = eye(d) - A * A' - B * B';
        Gt=C + A * (rotated_returns' * rotated_returns) * A' + B * Gt_prev * B';
    elseif strcmp(model, 'RDCC')
        % Additional computations for RDCC if necessary  
        % CAMBIAR PORQUE NO HICE RDCC TODAV�?A
        C = eye(d) - A * A' - B * B';
        Gt=C + A * (rotated_returns' * rotated_returns) * A' + B * Gt_prev * B';
        
    end
end 