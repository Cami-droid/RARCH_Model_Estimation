models = {'RBEKK', 'OGARCH', 'GOGARCH', 'RDCC'};
specifications = {'Scalar', 'Diagonal', 'CP'}; 
for i=1:4
    for j=1:3
     et=outputs(i,j).rotated_returns(end)
     Gt=outputs(i,j).Gt(:,:,end)
    end
end