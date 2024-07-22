function ll = ll_type(model, specification, outputs, t)
    thetaS = outputs.H_bar;    
    Dt = outputs.Dt;
    Ct = outputs.Ct;
    et = outputs.rotated_returns;
    Gt = outputs.Gt;
    d = outputs.d;

    fprintf('Calculando ll_type para el modelo %s en t=%d\n', model, t);
    fprintf('Gt size: %s\n', mat2str(size(Gt(:,:,t))));
    fprintf('et size: %s\n', mat2str(size(et)));
    disp('Gt:'); disp(Gt(:,:,t));

    delta = 1;
    U(delta) = 1;

    switch model
        case 'RBEKK'
            ll =  -0.5 * (d * log(2*pi) + log(det(Gt(:,:,t))) + et(t, :) * inv(Gt(:,:,t)) * et(t, :)');

        case 'OGARCH'
            ll =  -0.5 * (d * log(2*pi) + log(det(Gt(:,:,t))) + et(t, :) * inv(Gt(:,:,t)) * et(t, :)');

        case 'GOGARCH'
            ll =  -0.5 * (d * log(2*pi) + log(det(Gt(:,:,t))) + et(t, :) * U(delta)' * inv(Gt(:,:,t)) * U(delta) * et(t, :)');

        case 'RDCC'
            fprintf('Calculando RDCC, t=%d\n', t);
            fprintf('Dt size: %s\n', mat2str(size(Dt(:,:,t))));
            fprintf('Ct size: %s\n', mat2str(size(Ct(:,:,t))));

            ll_1 = -0.5 * (d * log(2*pi) + 2*log(det(Dt(:,:,t)))) + et(t, :) * inv(Dt(:,:,t)^2) * et(t, :)';
            ll_2 = -0.5 * (-et(t, :) * et(t, :)' + log(det(Ct(:,:,t))) + et(t, :) * inv(Ct(:,:,t)) * et(t, :)');

            ll = ll_1 + ll_2;
    end

    fprintf('Log-likelihood calculada: %f\n', ll);
end

