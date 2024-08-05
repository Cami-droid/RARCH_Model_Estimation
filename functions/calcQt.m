function [Qt, Qt_star , Ct] = calcQt(model,specification, outputs, thetaD)

    %  Initialization variables
    [i, j] = models_index(model, specification);
    rotated_returns = outputs.rotated_returns;
    T = outputs.T;
    d = outputs.d;
    thetaS=outputs.PI_bar;
    P_C = outputs.P;
    Gt =outputs.Gt;
    Lambda_C= outputs.Lambda_C;
    Dt = outputs.Dt;
    Id=eye(d);

    % Prealocate Qt and Ct for efficiency

    % Initialize Qt y Qt_star if i==3
    if i == 4 % 4 is the model index for RDCC
        % Calculate Qt y Ct for t = 0 
        initial_Qt_star=Id;
        Qt_star = calc_all_Gts(model,specification,outputs, thetaD, initial_Qt_star); %%% I use Gt because it's the same formula
        Qt = zeros(d, d, T);
        initial_Qt = P_C * sqrt(Lambda_C) * P_C' * initial_Qt_star * P_C * sqrt(Lambda_C) * P_C';
        
        Ct = zeros(d, d, T);
        initial_Ct = sqrt(initial_Qt .* Id) * initial_Qt * sqrt(inv(initial_Qt .* Id));

        % Calculate Qt and Ct for t >= 1
        for     t=1:T
        Qt(:, :, t) = P_C * sqrt(Lambda_C) * P_C' * Qt_star(:, :, t) * P_C * sqrt(Lambda_C) * P_C';
        Ct(:, :, t) = sqrt(inv(Qt(:, :, t) .* Id)) * Qt(:, :, t) * sqrt(inv(Qt(:, :, t) .* Id));
        end
    else
        Qt = [];
        Qt_star = [];
        Ct=[];
    end