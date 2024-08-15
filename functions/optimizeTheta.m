function [theta_opt, LL_total, Gt, VCV, Scores] = optimizeTheta(model, specification, outputs, initial_thetaD)
    % Function to optimize theta parameters for different models
    % Inputs:
    %   model          - The model name ('RBEKK', 'OGARCH', 'GOGARCH', 'RDCC')
    %   specification  - The specification type
    %   outputs        - The structure containing necessary data and parameters
    %   initial_thetaD - Initial thetaD parameters
    % Outputs:
    %   theta_opt      - Optimized theta parameters vector
    %   LL_total           - Final log-likelihood value
    %   Gt             - Estimated Gt matrices
    %   VCV            - Variance-covariance matrix
    %   Scores         - Scores

    % Extract necessary data
    d = outputs.d;
    et = outputs.rotated_returns;
    rt = outputs.returns;

    %  optimization options
    options = optimset('fmincon');
    options.Display = 'off'; % 'iter'
    options.Diagnostics = 'off'; % 'on'
    options.Algorithm ='interior-point';%'sqp'
    options.MaxIterations=100; 

    % Initialize outputs
    
    theta_opt = [];
    LL_total = [];
    Gt = [];
    VCV = [];
    Scores = [];

    % Switch case to call the appropriate optimization function based on the model
    try
        switch model
            case 'RBEKK'
                [theta_opt, LL_total, Gt, VCV, Scores] = rarch(et, 1, 1, specification, '2-stage',[],options);
                theta_opt=theta_opt'; %the rarch function parameter output is a column vector, we need a row vector 
                %rarch(data,p,q,type,method,startingVals,options)

            case 'OGARCH' 
            
                %   OUTPUTS:
                %   PARAMETERS   - Estimated parameters in the order:
                %                    OGARCH:
                %                    [vol(1) ... vol(K)]
                %                    GOGARCH:
                %                    [phi(1) ... phi(K(K-1)/2) vol(1) ... vol(K)]
                %                    where vol(i) = [alpha(i,1) ... alpha(i,P(i)) beta(i,1) ... beta(i,Q(i))]
                %                       as vol(i) are pairs alpha and beta, and we need first a vector of alphas and second a vector of betas
                %                       it's need to rearrange the output parameter vector 

                [theta_opt, LL_total, Gt, VCV, Scores] = gogarch_spec(et, 1, 1, [], model,[],options,specification);

                theta_alphas=theta_opt(1:2:end);
                theta_betas=theta_opt(2:2:end);
                theta_opt=[theta_alphas,theta_betas];
            case 'GOGARCH'

                [theta_opt, LL_total, Gt, VCV, Scores] = gogarch_spec(et, 1, 1, [], model,[],options,specification);

                %rearranging theta_opt

                theta_opt_phipart=theta_opt(1:(d*(d-1)/2));
                theta_opt_volpart=theta_opt(((d*d-d+2)/2):end);
                theta_alphas=theta_opt_volpart(1:2:end);
                theta_betas=theta_opt_volpart(2:2:end);
                theta_opt=[theta_opt_phipart,theta_alphas,theta_betas];

            
            case 'RDCC'
            [theta_opt, LL_total ,Gt, VCV, scores, ~]=rdcc(rt,[],1,1,1,1,2,'3-stage','None',[],options,specification);
            otherwise
                error('Unknown model: %s', model);
        end
    catch ME
        fprintf('Error optimizing model %s with specification %s: %s\n', model, specification, ME.message);
    end

    % Additional checks or calculations can be added here if needed

end 