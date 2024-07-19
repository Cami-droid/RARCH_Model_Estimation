function cond_cov = conditional_covariance(params, model)
    switch model
        case 'RBEKK'
            % Cálculo de covarianza condicional para RBEKK
        case 'OGARCH'
            % Cálculo de covarianza condicional para OGARCH
        case 'GOGARCH'
            % Cálculo de covarianza condicional para GOGARCH
        case 'RDCC'
            % Cálculo de covarianza condicional para RDCC
    end
    
    cond_cov = calculated_covariance;
end
