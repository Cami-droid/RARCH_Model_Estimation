function [A,B,C]=ABC(model, specification, outputs)
d=outputs.d
[i,j] = models_index(model, specification);


switch specification
case 'Scalar'
A=thetaD_opt(1)*eye(d)
B=thetaD_opt(2)*eye(d)
C=C = eye(d) - A * A' - B * B'
case 'Diagonal'
thetaD_alphas=thetaD_opt(1:d)
thetaD_betas=thetaD_opt(d+1:2*d)
A= diag(thetaD_alphas)
B= diag(thetaD_betas)
C=C = eye(d) - A * A' - B * B'

case 'CP'
lambda_cp=thetaD_opt(end)
thetaD_alphas=thetaD_opt(1:d)
thetaD_betas=sqrt(lambda_cp-thetaD_alphas)
A= diag(thetaD_alphas)
B= diag(thetaD_betas)
C = eye(d) - A * A' - B * B'
end
