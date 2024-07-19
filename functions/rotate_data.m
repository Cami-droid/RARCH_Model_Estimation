function [rotated_data, H_bar, Lambda, P] = rotate_data(returns, model)
    H_bar = cov(returns);

    disp('Unconditional Covariance Matrix H_bar');
    disp(H_bar)

    % Descomposicion en valores propios
    [P, Lambda] = eig(H_bar);

    % Matriz de eigenvectores P y matriz diagonal de eigenvalores Lambda
    disp('Eigenvectors matrix P:');
    disp(P);
    disp('Eigenvalues Diagnonal matrix Lambda:');
    disp(Lambda);

    switch model
        case {'RBEKK', 'GOGARCH', 'RDCC'}
            % Rotación para RBEKK, GOGARCH y RDCC
            H_bar_inv_sqrt = P * sqrt(inv(Lambda)) * P';

            e_t = returns*H_bar_inv_sqrt;

             
        case {'OGARCH'}
            % Rotación específica para OGARCH

            H_bar_inv_sqrt = P * sqrt(inv(Lambda)) * P';

            e_t = returns*sqrt(inv(Lambda)) * P';
       
    end
    
    % disp('Rotated returns e_t:');
    % disp(e_t);
    rotated_data = e_t;
end
