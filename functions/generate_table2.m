function generate_table(results, outputs)
    % Inicializar una celda para almacenar los datos de la tabla
    table_data = {};
    
    % Encabezados de la tabla
    headers = {'Model', 'Specification', 'ThetaD_opt', 'Fval', 'Qt', 'Qt_star'};
    
    % Agregar encabezados a la tabla
    table_data = [table_data; headers];
    
    % Iterar sobre los modelos y especificaciones
    [num_models, num_specs] = size(results);
    for i = 1:num_models
        for j = 1:num_specs
            % Extraer la información relevante de results y outputs
            model = results(i, j).model;
            specification = results(i, j).specification;
            thetaD_opt = mat2str(results(i, j).thetaD_opt);
            fval = results(i, j).fval;
            
            % Qt y Qt_star se almacenan como matrices 3D, se pueden promediar o sumar para un resumen
            Qt_mean = mean(outputs(i, j).Qt, 3);
            Qt_star_mean = mean(outputs(i, j).Qt_star, 3);
            
            % Convertir matrices a strings para ser agregadas a la tabla
            Qt_str = mat2str(Qt_mean);
            Qt_star_str = mat2str(Qt_star_mean);
            
            % Agregar los datos a la tabla
            table_data = [table_data; {model, specification, thetaD_opt, fval, Qt_str, Qt_star_str}];
        end
    end
    
    % Convertir la celda a tabla
    T = cell2table(table_data(2:end,:), 'VariableNames', table_data(1,:));
    
    % Especificar el nombre del archivo Excel
    filename = 'RARCH_Model_Results.xlsx';
    
    % Escribir la tabla en un archivo Excel
    writetable(T, filename);
    
    % Mostrar un mensaje de éxito
    fprintf('Tabla generada y guardada en %s\n', filename);
end
