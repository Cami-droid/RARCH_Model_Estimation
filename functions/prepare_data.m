function returns = prepare_data(log_returns,model)
    
    % Calcular la media de los log returns
    mean_log_returns = mean(log_returns);

    % Calcular los residuales (media cero)
    demeaned_returns = log_returns - mean_log_returns;
    
    % Preparar los datos según el modelo
    switch model
        case 'RBEKK'
            % Preparación específica para RBEKK.
            prepared_data=demeaned_returns

            
        case 'OGARCH'
            % Preparación específica para OGARCH

            prepared_data=demeaned_returns
        case 'GOGARCH'
            % Preparación específica para GOGARCH

            prepared_data=demeaned_returns
        case 'RDCC'
            % Preparación específica para RDCC

            prepared_data=standardize(log_returns)
    end
    
    % Devolver los datos preparados
    returns = prepared_data;
    
end
