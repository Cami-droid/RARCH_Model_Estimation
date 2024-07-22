function ll=ll_type(model,d,thetaS,Gt,et,i,Dt,Ct)
% ll_type(model,d,Gt,et,i,D,C)
% Dt and Ct have to be fulfilled only when computing RDCC model logLikelihood
delta=1;
U(delta)=1;

    switch model
        case 'RBEKK'

            ll =  -0.5 * (d * log(2*pi) + log(det(Gt)) + et(i-1, :) * inv(Gt) * et(i-1, :)');

        case 'OGARCH'

            ll =  -0.5 * (d * log(2*pi) + log(det(Gt)) + et(i-1, :) * inv(Gt) * et(i-1, :)');

        case 'GOGARCH'

            ll =  -0.5 * (d * log(2*pi) + log(det(Gt)) + et(i-1, :) * U(delta)'*inv(Gt) * U(delta)*et(i-1, :)');

        case 'RDCC'
            
            calcGt(thetaS, thetaD, rotated_returns, initial_Gt, model, specification)
            
            fprintf('estoy en la función ll_type me da error det(Dt), este es el size')
            size(Dt(:,:,i))
            Dt(:,:,i)
            ll_1= -0.5 * (d * log(2*pi) + 2*log(det(Dt(:,:,i))))+et(i-1, :)*Dt(:,:,i)^(-2)*et(i-1, :)' ;   %cambiar por parï¿½metro dependiente de t D(:,:,t)^(-2)
            ll_2= -0.5 * (-et(i-1, :)*et(i-1, :)' +log(det(Ct(:,:,i))) + et(i-1, :)*inv(Ct(:,:,i))*et(i-1, :)');

            ll=ll_1+ll_2;
    end

end