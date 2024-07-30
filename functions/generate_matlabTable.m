% Assuming that 'results' is the structure with a size of 4x3

models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};

% Define all possible parameter names without special characters
all_params = {'alpha', 'beta', 'alpha_11', 'alpha_22', 'beta_11', 'beta_22', 'lambda_cp', 'delta'};

% Initialize a cell array to store the table
thetaD_table = cell(numel(all_params) + 1, numel(models) * numel(specifications));

% Fill the table with transposed thetaD_opt vectors
for i = 1:4
    for j = 1:3
        % Get the vector of transposed parameters
        thetaD_opt = results(i, j).thetaD_opt';
        
        % Identify which parameters should be included for this
        % combination of model and specification
        param_idx = false(size(all_params));
        switch specifications{j}
            case 'Scalar'
                param_idx(1:2) = true; % alpha and beta
            case 'Diagonal'
                param_idx(3:6) = true; % alpha_11, alpha_22, beta_11, beta_22
            case 'CP'
                param_idx([1 2 7]) = true; % alpha, beta, lambda_cp
        end
        
        if i == 3
            param_idx(8) = true; % delta
        end
        
        % Fill the cells with the parameter values
        thetaD_table(param_idx, (i-1)*3 + j) = num2cell(thetaD_opt);
        
        % Add the likelihood value
        thetaD_table{numel(all_params) + 1, (i-1)*3 + j} = results(i, j).fval;
    end
end

% Create column names
col_names = cell(1, numel(models) * numel(specifications));
for i = 1:numel(models)
    for j = 1:numel(specifications)
        col_names{(i-1)*3 + j} = [models{i} '_' specifications{j}];
    end
end

% Add the row name for likelihood
all_params{end + 1} = 'likelihood';

% Create a table in MATLAB to visualize the results
thetaD_table_matlab = cell2table(thetaD_table, 'VariableNames', col_names, 'RowNames', all_params);

% Display the table
disp(thetaD_table_matlab);

