time_span=10
models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'}; 
for i=1:4
    for j=1:3
     et=outputs(i,j).rotated_returns(end)
     u=mean(outputs.i,j.returns)
     G_prev=outputs(i,j).Gt(:,:,end)
     G=calcG
     f_et= % f stands by forecast
    end
end

