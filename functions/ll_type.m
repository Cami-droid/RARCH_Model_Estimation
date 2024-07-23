function [ll results]= ll_type(model, specification, outputs, t)

    [i,j]=models_index(model,specification);
    thetaS = outputs.H_bar;    
    Dt = outputs.Dt;
    Ct = outputs.Ct;
    et = outputs.rotated_returns;
    rt=  outputs.returns;
    Gt = outputs.Gt;
    d = outputs.d;

    %fprintf('Calculando ll_type para el modelo %s en %d\n', model, t);
    
    %fprintf('Gt size: %s\n', mat2str(size(Gt(:,:,t))));
    %fprintf('et size: %s\n', mat2str(size(et)));
    %disp('Gt:'); disp(Gt(:,:,t));

    
    % Definir delta (puedes ajustar este valor según tus necesidades)
    delta = 1;
    
    % Obtener la matriz U(delta)
    U = U_delta(delta, d);
    
    if any(isnan(Gt(:,:,t)))
        error('Gt contiene valores NaN en t=%d', t);
    end

    switch model
        case 'RBEKK'
            ll =  -0.5 * (d * log(2*pi) + log(det(Gt(:,:,t))) + et(t, :) * inv(Gt(:,:,t)) * et(t, :)');

        case 'OGARCH'
            ll = -0.5 * (d * log(2 * pi) + log(det(Gt(:,:,t))) + et(t, :) * U * inv(Gt(:,:,t)) * U' * et(t, :)');

        case 'GOGARCH'
            ll = -0.5 * (d * log(2 * pi) + log(det(Gt(:,:,t))) + et(t, :) * U * inv(Gt(:,:,t)) * U' * et(t, :)');
        case 'RDCC'
            fprintf('Calculando RDCC, t=%d\n', t);
            %fprintf('Dt size: %s\n', mat2str(size(Dt(:,:,t))));
            %fprintf('Ct size: %s\n', mat2str(size(Ct(:,:,t))));

            ll_1 = -0.5 * (d * log(2*pi) + 2*log(det(Dt(:,:,t)))) + rt(t, :) * inv(Dt(:,:,t)^2) * rt(t, :)';
            ll_2 = -0.5 * (-et(t, :)*et(t, :)' + log(det(Ct(:,:,t))) + et(t, :)* inv(Ct(:,:,t))*et(t, :)');

            ll = ll_1 + ll_2;
    end

    
end

function U = U_delta(delta, d)
    % U_delta returns a diagonal matrix U with diagonal entries equal to delta
    U = delta * eye(d);
end
