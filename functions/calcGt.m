function Gt = calcGt(model, specification, outputs, thetaD, Gt_prev,t)
    if t>=2
        et = outputs.rotated_returns(t-1,:);
    else
        et=0; % backCast
    end
    thetaS = outputs.H_bar;
    d=outputs.d;

    if strcmp(model, 'GOGARCH')
        delta = thetaD(end); % when the model is GOGARCH δ is the last element of thetaD vector;
    else
        delta=1;  % put 1 to don't affect   
    end

    % Extract dynamic parameters based on specification
    if strcmp(specification, 'Scalar')

        alpha_vec = thetaD(1) * ones(1, d); 
        beta_vec = thetaD(2) * ones(1, d);

    elseif strcmp(specification, 'Diagonal')

        alpha_vec = thetaD(1:d);
        beta_vec = thetaD(d+1:2*d);

    elseif strcmp(specification, 'CP')

        alpha_vec = thetaD(1:d) ;
        lambda_cp = thetaD(d+1) ;
        beta_vec =(lambda_cp*ones(1,d) - alpha_vec); %% thetaD(d+1) is lambda and the thetaD(1:d) are alphas, all term root squared

    end
    
    A = diag(sqrt(alpha_vec));
    B = diag(sqrt(beta_vec));
    
    % Small regularization term to avoid singularity

    %reg_term = 1e-6 * eye(d); % regularization term
    C = eye(d) - A * A' - B * B';

    % Compute Gt based on the model

    if strcmp(model, 'RBEKK')

        
        Gt = C + (A * A').*(et' * et)  + (B * B').* Gt_prev; % Adding regularization term;
      
    elseif  strcmp(model, 'OGARCH')

        
        Gt = C + A * A'.*(et' * et)  + B * B'.* Gt_prev; % Adding regularization term

    elseif  strcmp(model, 'GOGARCH')
        et=delta*et;

               
        %%%%%%%%%%%%%%% my formulas generate this algebraic expression, reminder:check why ll goes too high%%%

        Gt = C + A  * A'.*(et' * et)  + B * B'.* Gt_prev; % Adding regularization term

    elseif strcmp(model, 'RDCC')

        Gt = C +A * A'.*(et' * et)  +B * B'.* Gt_prev ;% Adding regularization term
    end

    % Ensure Gt is symmetric
    Gt = (Gt + Gt') / 2;
end
