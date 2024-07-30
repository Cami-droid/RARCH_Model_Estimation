% Definir modelos y especificaciones
models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};

% Definir nombres de parámetros
marginal_params = {'alpha', 'beta'};
dynamic_params = {'alpha', 'beta', 'alpha_11', 'alpha_22', 'beta_11', 'beta_22', 'lambda_cp', 'delta'};
d = 2; %%%%%%%%%%% probar sacandolo%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ll_decomposition_params = [arrayfun(@(x) ['LL_', num2str(x)], 1:d, 'UniformOutput', false), {'Copula_LL', 'Total_LL'}];

% Calcular número de filas por sección
num_marginal_params = numel(marginal_params);
num_dynamic_params = numel(dynamic_params);
num_ll_decomposition_params = numel(ll_decomposition_params);

% Inicializar matriz de celdas para la tabla completa
total_rows = num_marginal_params + num_dynamic_params + num_ll_decomposition_params + 3;
total_cols = numel(models) * numel(specifications);
complete_table = cell(total_rows, total_cols + 1);

% Rellenar títulos de subtabla
complete_table{1, 1} = 'Marginal Parameters';
complete_table{num_marginal_params + 2, 1} = 'Dynamic Parameters';
complete_table{num_marginal_params + num_dynamic_params + 3, 1} = 'LL Decomposition';

% Rellenar encabezados de fila para "Marginal Parameters"
for i = 1:num_marginal_params
    complete_table{i + 1, 1} = marginal_params{i};
end

% Rellenar encabezados de fila para "Dynamic Parameters"
for i = 1:num_dynamic_params
    complete_table{i + num_marginal_params + 2, 1} = dynamic_params{i};
end

% Rellenar encabezados de fila para "LL Decomposition"
for i = 1:num_ll_decomposition_params
    complete_table{i + num_marginal_params + num_dynamic_params + 3, 1} = ll_decomposition_params{i};
end

% Rellenar datos de results en la tabla
for i = 1:4
    for j = 1:3
        thetaD_opt = results(i, j).thetaD_opt';
        
        % "Marginal Parameters" (suponiendo que los primeros dos parámetros son marginales)
        complete_table(2:3, (i-1)*3 + j + 1) = num2cell(thetaD_opt(1:2));
        
        % "Dynamic Parameters"
        param_idx = false(size(dynamic_params));
        switch specifications{j}
            case 'Scalar'
                param_idx(1:2) = true;
            case 'Diagonal'
                param_idx(3:6) = true;
            case 'CP'
                param_idx([1 2 7]) = true;
        end
        if i == 3
            param_idx(8) = true;
        end
        
        complete_table(find(param_idx) + num_marginal_params + 2, (i-1)*3 + j + 1) = num2cell(thetaD_opt);
        
        % "LL Decomposition"
        ll_values = [NaN(d, 1); NaN; results(i, j).fval];
        complete_table(num_marginal_params + num_dynamic_params + 4:end, (i-1)*3 + j + 1) = num2cell(ll_values);
    end
end

% Redondear valores numéricos a 3 decimales
num_decimals=3;
for i = 2:size(complete_table, 1)
    for j = 2:size(complete_table, 2)
        if ~isempty(complete_table{i, j}) && isnumeric(complete_table{i, j})
            complete_table{i, j} = round(complete_table{i, j}, num_decimals);
        end
    end
end

% Crear nombres únicos de columna para la tabla completa
col_names = cell(1, total_cols + 1);
col_names{1} = 'Parameter';
for i = 1:numel(models)
    for j = 1:numel(specifications)
        col_names{(i-1)*3 + j + 1} = [models{i}, '_', specifications{j}];
    end
end

% Convert the complete_table to a MATLAB table
complete_table_matlab = cell2table(complete_table, 'VariableNames', col_names);

% Display the complete table
disp(complete_table_matlab);


% Define the path and filename for the Excel file
results_dir = 'D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\results';
excel_file = fullfile(results_dir, 'complete_table.xlsx');

% Write the complete table to the Excel file
writetable(complete_table_matlab, excel_file, 'WriteVariableNames', true, 'WriteRowNames', false);

% Display a success message
disp(['Table has been successfully exported to ', excel_file]);
