% Main file to run all the functions and generate the final table

% Add 'functions' folder to the path
addpath('functions'); 

models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};

I = length(models);
J = length(specifications);

% Initialize results and outputs structures
results(I, J) = struct('model', [], 'specification', [], 'thetaM', [], 'theta_vec', [], 'fval', [], 'Qt', [], 'Qt_star', [], 'L', [], 'LL_marginal', [], 'LL_copula', []);
outputs(I, J) = struct('model', [], 'specification', [], 'P', [], 'Lambda', [], 'H_bar', [], 'Gt', [], 'Gt_artesanal', [], 'VCV', [], 'Scores', [], 'returns', [], 'initials_thetaD', [], 'rotated_returns', [], 'Dt', [], 'Ct', [], 'I', [], 'J', [], 'd', [], 'T', [], 'L', [], 'LL_marginal', [], 'LL_copula', []);

alpha_init = 0.02;
beta_init = 0.85;

initials_thetaD = { [alpha_init, beta_init], [alpha_init*ones(1,d), beta_init*ones(1,d)], [alpha_init*ones(1,d), alpha_init + beta_init] };
initial_delta = 0.1; % Only for 'GOGARCH' models

% Prepare the outputs structure
for i = 1:I
    for j = 1:J
        outputs(i,j).I = I;
        outputs(i,j).J = J;
        outputs(i,j).d = d;
        outputs(i,j).T = T;
        outputs(i,j).Gt = zeros(d, d, T); % T+1 because the first matrix is index 0 in theory
        outputs(i,j).initial_Gt = eye(d);
        outputs(i,j).initials_thetaD = initials_thetaD{j};

        if strcmp(models{i}, 'GOGARCH')
            outputs(i,j).initials_thetaD = [initials_thetaD{j}, initial_delta];
        end
    end
end

% Estimation and Optimization
for i = 1:I
    model = models{i};
    fprintf('*********************************** Estimating model: %s ****************************************\n', model);
    
    for j = 1:J
        specification = specifications{j};
        outputs(i,j).model = model;
        outputs(i,j).specification = specification;

        % Load and prepare data
        [outputs(i,j).returns, outputs(i,j).Dt, results(i,j).thetaM] = prepare_data(model, outputs(i,j), log_returns);

        % Rotate data
        [outputs(i,j).rotated_returns, outputs(i,j).H_bar, outputs(i,j).Lambda, outputs(i,j).P] = rotate_data(outputs(i,j), model);

        fprintf('Estimating model: %s with specification: %s\n', model, specification);
        fprintf('Initial parameters (thetaD): %s\n', mat2str(outputs(i,j).initials_thetaD));

        % Optimize logLikelihood
        initial_thetaD = outputs(i,j).initials_thetaD;
        [results(i,j).theta_vec, results(i,j).fval, outputs(i,j).Gt, outputs(i,j).VCV, outputs(i,j).Scores] = optimizeTheta(model, specification, outputs(i,j), initial_thetaD);

        fprintf('Optimal parameters (theta): %s\n', mat2str(results(i,j).theta_vec));
        fprintf('LogLikelihood value: %f\n', results(i,j).fval);

        % Calculate Qt, Qt_star and Ct
        [outputs(i,j).Qt, outputs(i,j).Qt_star, outputs(i,j).Ct] = calcQt(model, specification, outputs(i,j), results(i,j).theta_vec);

        % Store the results
        results(i,j).model = model;
        results(i,j).specification = specification;
    end
end

% Generate the final table
generate_matlabTable;
