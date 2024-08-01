function [theta_vec, fval, Gt, VCV,scores] = optimizeThetaD(model, specification, outputs, thetaD_initial)
%%% fval: log-verosimilitud final
    d=outputs.d;
    et=outputs.rotated_returns;
    rt=outputs.returns;
    

    switch model
        case 'RBEKK'
        
        [theta_vec,fval,Gt,VCV,Scores] = rarch(et,1,1,'Scalar','2-stage');

        [theta_vec,fval,Gt,VCV,Scores] = rarch(et,1,1,'Diagonal','2-stage');

        [theta_vec,fval, Gt,VCV,Scores] = rarch(et,1,1,'CP','2-stage');

        case 'OGARCH' || 'GOGARCH'

        [theta_vec, fval, Gt, VCV,Scores] = gogarch(e_t, p, q,[], model);

        
        case 'RDCC'

        [theta_vec, fval, Gt, VCV,Scores] = gogarch(e_t, p, q,[], model);
    end




