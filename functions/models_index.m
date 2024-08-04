function [i, j] = models_index(model, specification)
    % Define model names and specificacions
    models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
    specifications = {'Scalar', 'Diagonal', 'CP'};
    
    % To find the model index
    i = find(strcmp(models, model));
    
    % to fin the specification index
    j = find(strcmp(specifications, specification));
    
    % Verify both indexes are valid
    if isempty(i)
        error('Invalid model.');
    end
    
    if isempty(j)
        error('Invalid specification.');
    end
end
