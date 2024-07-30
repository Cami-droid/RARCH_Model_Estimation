function [rotated_data, H_bar, Lambda, P] = rotate_data(outputs, model)
    T=outputs.T;
    returns=outputs.returns;
            
    H_bar = (1/T).*returns'*returns;
    % Descomposicion en valores propios
    [P, Lambda] = eig(H_bar);
    H_bar_inv_sqrt = P * sqrt(inv(Lambda)) * P';
    
    switch model
        case {'RBEKK'}

            % Rotación para RBEKK y GOGARCH
            e_t = returns*H_bar_inv_sqrt;

             
        case {'OGARCH'}

            % Specific rotation for OGARCH;
            e_t = returns*H_bar_inv_sqrt * P';

        case {'GOGARCH'}
            
            % Specific rotation for RBEKK and GOGARCH
            e_t = returns*H_bar_inv_sqrt;
            
        case {'RDCC'}
            PI_bar=H_bar; %% porque los returns que vienen del prepare_data están estandarizados ya para RDCC
            %%% solo estoy dandole los nombres que están en el paper para evitar confusiones
            
            %disp('Unconditional Covariance Matrix PI_bar valid for RDCC');
            %disp(PI_bar);
            % Matriz de eigenvectores P y matriz diagonal de eigenvalores Lambda
            %disp('Eigenvectors matrix P:');
            %disp(P);
            %disp('Eigenvalues Diagnonal matrix Lambda:');
            %disp(Lambda);
            
            % specific rotation for RDCC
            PI_bar_inv_sqrt = P * sqrt(inv(Lambda)) * P';        
            e_t = returns*PI_bar_inv_sqrt;
            
       
    end
    % disp('Rotated returns e_t:');
    % disp(e_t);
    rotated_data = e_t;
end
