% Define model an specifications
models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};

% initialize parameters

% Define parameter names

marginal_params = {};
for i = 1:d
    marginal_params{end+1} = sprintf('alpha_%d', i);
    marginal_params{end+1} = sprintf('beta_%d', i);
end

%%  dynamic_params: dynamic parameter names that change according to d


dynamic_params = {'alpha', 'beta'}; % initializing cell with scalar dynamic parameters 

% Add the  parameters alpha and beta for each dimension
for i = 1:d
    dynamic_params{end+1} = sprintf('alpha_%d', i);
end
for i = 1:d
    dynamic_params{end+1} = sprintf('beta_%d', i);
end

% Add the parameters  lambda_cp and delta
dynamic_params{end+1} = 'lambda_cp';
dynamic_params{end+1} = 'delta';

% Show results
disp(dynamic_params);


%% 

ll_decomposition_params = [arrayfun(@(x) ['LL_', num2str(x)], 1:d, 'UniformOutput', false), {'Copula_LL', 'Total_LL'}];

% Calculate number of rows in each section
num_marginal_params = numel(marginal_params);
num_dynamic_params = numel(dynamic_params);
num_ll_decomposition_params = numel(ll_decomposition_params);

% Inicialize a cell matrix for the complete table
total_rows = num_marginal_params + num_dynamic_params + num_ll_decomposition_params + 3; % 3 subtable titles
total_cols = numel(models) * numel(specifications);
complete_table = cell(total_rows, total_cols + 1);

% fill subtable titles 
complete_table{1, 1} = 'Marginal Parameters';
complete_table{num_marginal_params  + 2, 1} = 'Dynamic Parameters';
complete_table{num_marginal_params + num_dynamic_params + 3, 1} = 'LL Decomposition';

% fill row headers for  "Marginal Parameters"

for i = 1:num_marginal_params
    complete_table{i + 1, 1} = [marginal_params{i}];
end


% fill row headers for "Dynamic Parameters"
for i = 1:num_dynamic_params
    complete_table{i + num_marginal_params + 2, 1} = dynamic_params{i};
end

% fill row headers for "LL Decomposition"
for i = 1:num_ll_decomposition_params
    complete_table{i + num_marginal_params + num_dynamic_params + 3, 1} = ll_decomposition_params{i};
end

% fill results data in the table
for i = 1:4
    for j = 1:3
        thetaD = results(i, j).thetaD;
        thetaM = results(i, j).thetaM; % Columns vector [alpha1;beta1;alpha2;beta2;...;alphad;betad]

        % "Marginal Parameters"
        if ~isempty(thetaM) && numel(thetaM) >= 2 * d
            for k = 1:d
                alpha_idx = 2 * (k - 1) + 1;
                beta_idx = 2 * (k - 1) + 2;
                complete_table{2 * (k - 1) + 2, (i - 1) * 3 + j + 1} = thetaM(alpha_idx);
                complete_table{2 * (k - 1) + 3, (i - 1) * 3 + j + 1} = thetaM(beta_idx);
            end
        end
        
        % "Dynamic Parameters"
        param_idx = false(size(dynamic_params)); % zero matriz with same size as dynamic params names
        switch specifications{j}
            case 'Scalar'
                % fills the two first rows 

                param_idx(1:2) = true; 
            case 'Diagonal'

                % fills starting from 3rd row
                param_idx(3:(2*d+2)) = true; %p=1 q=1 so 2=1+1
            case 'CP'
                % fills starting from 3rd row
                % the end -1 param_idx is the lambda_cp place

                param_idx([3:(d+2),end-1]) = true; 
        end
        if i == 3
             % the last param_idx is the delta place
            param_idx(end) = true;
        end

        % search the right place to write thetaD vector

        % number 2 stands for the 2 subtable title's places
        complete_table([num_marginal_params+find(param_idx) + 2], (i - 1) * 3 + j + 1) = num2cell(thetaD);
        
        % "LL Decomposition"
        ll_values = [results(i, j).LL_marginal; results(i, j).LL_copula; results(i, j).LL_total];
        complete_table((num_marginal_params  + num_dynamic_params + 4):end, (i - 1) * 3 + j + 1) = num2cell(ll_values);
    end
end

% round numerical values to 3 digits
num_decimals = 3;
for i = 2:size(complete_table, 1)
    for j = 2:size(complete_table, 2)
        if ~isempty(complete_table{i, j}) && isnumeric(complete_table{i, j})
            complete_table{i, j} = round(complete_table{i, j}, num_decimals);
        end
    end
end

% Create unique column numbers  for the complete table
col_names = cell(1, total_cols + 1);
col_names{1} = 'Parameter';
for i = 1:numel(models)
    for j = 1:numel(specifications)
        col_names{(i - 1) * 3 + j + 1} = [models{i}, '_', specifications{j}];
    end
end

% Convert complete_table into a matlab table
complete_table_matlab = cell2table(complete_table, 'VariableNames', col_names);

% show the complete table 
disp(complete_table_matlab);

% Define the path and file name for the Excel file
% Define the relative path to the results directory
results_dir = fullfile(pwd, 'results');
file_name = sprintf('complete_%s.xlsx', Task); % Incluye el valor de Task en el nombre del archivo
excel_file = fullfile(results_dir, file_name);




% Write the complete table in the Excel file
writetable(complete_table_matlab, excel_file, 'WriteVariableNames', true, 'WriteRowNames', false);

% Save structs 'outputs' and 'results' in MAT files

% Save 'outputs' variable in a file named 'outputs.mat' in the results folder
save(fullfile(results_dir, 'outputs.mat'), 'outputs');

% Save 'results' variable in a file named 'results.mat' in the results folder
save(fullfile(results_dir, 'results.mat'), 'results');


fprintf('Outputs and results saved successfully in the "results" folder.\n');

% Show success message
disp(['Table has been successfully exported to ', excel_file]);
