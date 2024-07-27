function thetaD_table_matlab = generate_pdfTable(results)
    models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
    param_types = {'Scalar', 'Diagonal', 'CP'};
    all_params = {'\alpha', '\beta', '\alpha_{11}', '\alpha_{22}', '\beta_{11}', '\beta_{22}', '\lambda_{CP}', '\delta'};

    % Inicializar la tabla con celdas vacías
    thetaD_table = cell(length(all_params), length(models) * length(param_types));
    
    for i = 1:length(models)
        for j = 1:length(param_types)
            idx = (i-1) * length(param_types) + j;
            thetaD_opt = results(i, j).thetaD_opt;
            switch param_types{j}
                case 'Scalar'
                    thetaD_table{1, idx} = thetaD_opt(1); % alpha
                    thetaD_table{2, idx} = thetaD_opt(2); % beta
                case 'Diagonal'
                    thetaD_table{3, idx} = thetaD_opt(1); % alpha_11
                    thetaD_table{4, idx} = thetaD_opt(2); % alpha_22
                    thetaD_table{5, idx} = thetaD_opt(3); % beta_11
                    thetaD_table{6, idx} = thetaD_opt(4); % beta_22
                case 'CP'
                    thetaD_table{1, idx} = thetaD_opt(1); % alpha
                    thetaD_table{2, idx} = thetaD_opt(2); % beta
                    thetaD_table{7, idx} = thetaD_opt(3); % lambda_CP
            end
        end
    end

    % Nombres de las columnas
    col_names = {};
    for i = 1:length(models)
        for j = 1:length(param_types)
            col_names{end+1} = sprintf('%s_%s', models{i}, param_types{j});
        end
    end

    % Crear la tabla en MATLAB
    thetaD_table_matlab = cell2table(thetaD_table, 'VariableNames', col_names, 'RowNames', all_params);

    % Guardar la tabla en un archivo .txt en el directorio results
    results_dir = 'D:\Documents\TRABAJO\Upwork\Rarch_model\work\RARCH_Model_Estimation\results';
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);knlkj
    end
    writetable(thetaD_table_matlab, fullfile(results_dir, 'thetaD_table.txt'), 'WriteRowNames', true, 'Delimiter', 'tab');
end
