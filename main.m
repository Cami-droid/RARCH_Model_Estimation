% Main file to run all the functions and generate the final table

% Add 'functions' folder to the path
addpath('functions'); 

models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};
             
I = length(models);
J = length(specifications);

results(I, J) = struct('model', [], 'specification', [], 'thetaM',[],'theta_vec', [], 'fval', [], 'Qt', [], 'Qt_star', [],'L',[],'LL_marginal',[],'LL_copula',[]);
outputs(I, J) = struct('model', [], 'specification', [], 'P', [], 'Lambda', [], 'H_bar', [], 'Gt', [],'Gt_artesanal',[],'VCV',[],'Scores',[],'returns', [], 'initials_thetaD',[],'rotated_returns', [], 'Dt', [], 'Ct', [], 'I', [], 'J', [], 'd', [], 'T', [],'L',[],'LL_marginal',[],'LL_copula',[]);
alpha_init=0.02;
beta_init=0.85;

initials_thetaD= { [alpha_init, beta_init], [alpha_init*ones(1,d),beta_init*ones(1,d)],[alpha_init*ones(1,d), alpha_init+beta_init]};

initial_delta=0.1; % only for 'GOGARCH' models
for i = 1:I
    
    for j = 1:J
        outputs(i,j).I = I;
        outputs(i,j).J = J;
        outputs(i,j).d = d;
        outputs(i,j).T = T;
        outputs(i,j).Gt = zeros(d, d, T + 1); % T+1 because the first matrix is index 0 in theory
        outputs(i,j).initial_Gt = eye(d);
        outputs(i,j).initials_thetaD=initials_thetaD{j};

        if i==3
        outputs(i,j).initials_thetaD=[initials_thetaD{j} initial_delta];
        end
                        
    end
end
clear delta

%% 
 
for i = 1:I
    model = models{i};
    fprintf('*********************************** Estimating model: %s ****************************************\n', model);
    for j = 1:J
        % Load and prepare data
        [outputs(i,j).returns, outputs(i,j).Dt, results(i,j).thetaM] = prepare_data(model, outputs(i,j), log_returns);
        
        % Rotate data
        [outputs(i,j).rotated_returns, outputs(i,j).H_bar, outputs(i,j).Lambda, outputs(i,j).P] = rotate_data(outputs(i,j), model);

         %   disp('Unconditional Covariance Matrix H_bar valid for RBEKK, OGARCH and GOGARCH');
         %   disp(outputs(i,j).H_bar);
         %   % Matriz de eigenvectores P y matriz diagonal de eigenvalores Lambda
         %   disp('Eigenvectors matrix P:');
         %   disp(outputs(i,j).P);
         %   disp('Eigenvalues Diagnonal matrix Lambda:');
         %   disp(outputs(i,j).Lambda);

    end

    for j = 1:J
        specification = specifications{j};
        outputs(i,j).model = model;
        outputs(i,j).specification = specification;

        fprintf('Estimating model: %s\n', model);
        fprintf('Parameter specification: %s\n', specification);
        
        % Optimize logLikelihood
        initial_thetaD = outputs(i,j).initials_thetaD;

        fprintf('The initials thetaDs are: %s\n', mat2str(outputs(i,j).initials_thetaD));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [results(i,j).theta_vec, results(i,j).fval, outputs(i,j).Gt, outputs(i,j).VCV, outputs(i,j).Scores] = optimizeTheta(model, specification, outputs(i,j), initial_thetaD)
        
        % Storing LL_marginals in results
% 
%         results(i,j).LL_marginal=sum(outputs(i,j).L_marginal)';
%         results(i,j).LL_copula=results(i,j).fval-sum(results(i,j).LL_marginal);
% 

        fprintf('The optimal vector theta found is: %s\n', mat2str(results(i,j).theta_vec));
        fprintf('LogLikelihood value: %s\n', mat2str(results(i,j).fval));

        % calculate Gt at the optimum theta_vec

        thetaD=results(i,j).theta_vec(((d+1)*d/2)+1:end)'
        
        % Calculate Qt, Qt_star and Ct
        
        [outputs(i,j).Qt, outputs(i,j).Qt_star, outputs(i,j).Ct] = calcQt(model, specification, outputs(i,j), results(i,j).theta_vec);
        
        clear thetaD
        % Store the results
        results(i,j).model = model;
        results(i,j).specification = specification;
    end
end
%% Table Generation
generate_matlabTable;


