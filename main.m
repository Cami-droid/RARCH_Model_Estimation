% Main file to run all the functions and generate the final table
clear;clc;

% Add 'functions' folder to the path
addpath('functions'); 

data = readtable('D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\data\stock_prices_28_KFT_UTX_1.csv');

% Extract data from relevant columns
dates = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
AA = data.AA; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMBIAR LO HAGO SOLO PARA REVISAR ERRORES
XOM = data.XOM;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMBIAR LO HAGO SOLO PARA REVISAR ERRORES
%%%%%%%%%%%%%%%%%%%%%REMEMBER THAT IT LOOSES ONE POSITION AFTER DIFFERENTIATION

% Calculate log returns
log_returns_AA = diff(log(AA)) * 100;
log_returns_XOM = diff(log(XOM)) * 100;

% Combine log returns in a matriz
log_returns = [log_returns_AA, log_returns_XOM];

d = size(log_returns, 2);
T = size(log_returns, 1);

% Assure that the date vector equalizes log returns' length
dates = dates(2:end);

models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};

             
I = length(models);
J = length(specifications);

results(I, J) = struct('model', [], 'specification', [], 'thetaD_opt', [], 'fval', [], 'Qt', [], 'Qt_star', [],'L',[]);
outputs(I, J) = struct('model', [], 'specification', [], 'P', [], 'Lambda', [], 'H_bar', [], 'Gt', [], 'Passenger_Gt',[],'returns', [], 'initials_thetaD',[],'rotated_returns', [], 'Dt', [], 'Ct', [], 'I', [], 'J', [], 'd', [], 'T', []);

    
for i = 1:I
    for j = 1:J
        outputs(i,j).I = I;
        outputs(i,j).J = J;
        outputs(i,j).d = d;
        outputs(i,j).T = T;
        outputs(i,j).Gt = zeros(d, d, T + 1); % T+1 because the first matrix is index 0 in theory
        outputs(i,j).initial_Gt = eye(d);
        outputs(i,j).initials_thetaD = { [0.01 0.3]      , [0.05 0.05 0.1 0.1]       , [0.05 0.05 0.2]};
        % Only for GOGARCH
        outputs(3,j).initial_delta=1;
    end
end
%% 

for i = 1:I
    model = models{i};
    fprintf('Estimating model: %s\n', model);
    for j = 1:J
        % Load and prepare data
        [outputs(i,j).returns, outputs(i,j).Dt] = prepare_data(model, outputs(i,j), log_returns);
        
        % Rotate data
        [outputs(i,j).rotated_returns, outputs(i,j).H_bar, outputs(i,j).Lambda, outputs(i,j).P] = rotate_data(outputs(i,j).returns, model);
    end

    for j = 1:J
        specification = specifications{j};
        outputs(i,j).model = model;
        outputs(i,j).specification = specification;

        fprintf('Estimating model: %s\n', model);
        fprintf('Parameter specification: %s\n', specification);
        
        % Optimize logLikelihood
        thetaD_initial = outputs(i,j).initials_thetaD{j};
        fprintf('The initial thetaD are: %s\n', mat2str(outputs(i,j).initials_thetaD{j}));
        
        % Estimate thetaD parameters
        outputs(i,j).Gt= calc_all_Gts(model, specification, outputs(i,j), outputs(i,j).initials_thetaD{j}, outputs(i,j).initial_Gt);
        
        [results(i,j).thetaD_opt, results(i,j).fval, exitflag, output] = optimizeThetaD(model, specification, outputs(i,j), thetaD_initial);
        
        disp(results(i,j).thetaD_opt)

        % calculate Gt at the optimum thetaD_opt

        Id=eye(d)
        output(i,j).Gt=calcGt(model, specification, outputs(i,j),results(i,j).thetaD_opt,Id)
        
        % Calculate Qt, Qt_star and Ct
        
        [outputs(i,j).Qt, outputs(i,j).Qt_star outputs(i,j).Ct] = calcQt(model, specification, outputs(i,j), results(i,j).thetaD_opt);
        
        % Store the results
        results(i,j).model = model;
        results(i,j).specification = specification;
    end
end

% Generate Table
generate_table(results);
