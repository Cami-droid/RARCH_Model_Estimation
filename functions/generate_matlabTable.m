%Assuming that 'results' is the structure that you have with a size of 4x3

models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'};

% Definir los nombres de todos los par�metros posibles sin caracteres especiales
all_params = {'alpha', 'beta', 'alpha_11', 'alpha_22', 'beta_11', 'beta_22', 'lambda_cp', 'delta'};

% Inicializar una celda para almacenar la tabla
thetaD_table = cell(numel(all_params), numel(models) * numel(specifications));

% Rellenar la tabla con los vectores thetaD_opt transpuestos
for i = 1:4
    for j = 1:3
        % Obtener el vector de par�metros transpuesto
        thetaD_opt = results(i, j).thetaD_opt';
        
        % Identificar qu� par�metros deben ser llenados para esta combinaci�n de modelo y especificaci�n
        param_idx = false(size(all_params));
        switch specifications{j}
            case 'Scalar'
                param_idx(1:2) = true; % alpha y beta
            case 'Diagonal'
                param_idx(3:6) = true; % alpha_11, alpha_22, beta_11, beta_22
            case 'CP'
                param_idx([1 2 7]) = true; % alpha, beta, lambda_cp
        end
        
        % Llenar las celdas correspondientes con los valores de los par�metros
        thetaD_table(param_idx, (i-1)*3 + j) = num2cell(thetaD_opt);
    end
end

% Crear nombres de las columnas
col_names = cell(1, numel(models) * numel(specifications));
for i = 1:numel(models)
    for j = 1:numel(specifications)
        col_names{(i-1)*3 + j} = [models{i} '_' specifications{j}];
    end
end

% Crear una tabla en MATLAB para visualizar los resultados
thetaD_table_matlab = cell2table(thetaD_table, 'VariableNames', col_names, 'RowNames', all_params);

% Mostrar la tabla
disp(thetaD_table_matlab);
