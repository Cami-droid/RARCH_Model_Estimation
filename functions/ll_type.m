function ll = ll_type(model, specification, outputs, t)
    t_count=t+1 % para computar en matlab y no confundirse con el verdadero t
    thetaS = outputs.H_bar;    
    Dt = outputs.Dt;
    Ct = outputs.Ct;
    et = outputs.rotated_returns;
    Gt = outputs.Gt;
    d = outputs.d;

    fprintf('Calculando ll_type para el modelo %s en t_count = t+1=%d\n', model, t);
    
    fprintf('Gt size: %s\n', mat2str(size(Gt(:,:,t_count))));
    fprintf('et size: %s\n', mat2str(size(et)));
    disp('Gt:'); disp(Gt(:,:,t));

    delta = 1;
    U(delta) = 1;
    

    if any(isnan(Gt(:,:,t)))
        error('Gt contiene valores NaN en t=%d', t);
    end

    switch model
        case 'RBEKK'
            ll =  -0.5 * (d * log(2*pi) + log(det(Gt(:,:,t_count))) + et(t, :) * inv(Gt(:,:,t_count)) * et(t, :)');

        case 'OGARCH'
            ll =  -0.5 * (d * log(2*pi) + log(det(Gt(:,:,t_count))) + et(t, :) * inv(Gt(:,:,t)) * et(t, :)');

        case 'GOGARCH'
            ll =  -0.5 * (d * log(2*pi) + log(det(Gt(:,:,t_count))) + et(t, :) * U(delta)' * inv(Gt(:,:,t_count)) * U(delta) * et(t, :)');

        case 'RDCC'
            fprintf('Calculando RDCC, t=%d\n', t);
            fprintf('Dt size: %s\n', mat2str(size(Dt(:,:,t_count))));
            fprintf('Ct size: %s\n', mat2str(size(Ct(:,:,t_count))));

            ll_1 = -0.5 * (d * log(2*pi) + 2*log(det(Dt(:,:,t_count)))) + et(t, :) * inv(Dt(:,:,t_count)^2) * et(t, :)';
            ll_2 = -0.5 * (-et(t, :) * et(t, :)' + log(det(Ct(:,:,t_count))) + et(t, :) * inv(Ct(:,:,t_count)) * et(t, :)');

            ll = ll_1 + ll_2;
    end

    fprintf('Log-likelihood calculada: %f\n', ll);
end
