function [theta_vec, fval, Gt, VCV,Scores] = optimizeTheta(model, specification, outputs, thetaD_initial)
%%% fval: log-verosimilitud final
    d=outputs.d;
    et=outputs.rotated_returns;
    rt=outputs.returns;
    

    switch model
        case 'RBEKK'
        
        [theta_vec,fval,Gt,VCV,Scores] = rarch(et,1,1,specification,'2-stage');

        case 'OGARCH' || 'GOGARCH'

        [theta_vec, fval, Gt, VCV,Scores] = gogarch(et, 1, 1,[], model);

        
        case 'RDCC'

        [theta_vec, fval, Gt, VCV,Scores] = dcc(et,[], 1, [],1,1,0,1,2,'3-stage','Scalar',[],[]);
        %[parameters, ll ,Ht, VCV, scores, diagnostics]=dcc(data,dataAsym,m,l,n,p,o,q,gjrType,method,composite,startingVals,options)
    end




