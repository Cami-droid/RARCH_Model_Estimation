function [i, j] = models_index(model, specification)
    % Definir los nombres de los modelos y especificaciones
    models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
    specifications = {'Scalar', 'Diagonal', 'CP'};
    
    % Encontrar el índice del modelo
    i = find(strcmp(models, model));
    
    % Encontrar el índice de la especificación
    j = find(strcmp(specifications, specification));
    
    % Verificar que ambos índices sean válidos
    if isempty(i)
        error('Modelo no válido.');
    end
    
    if isempty(j)
        error('Especificación no válida.');
    end
end
