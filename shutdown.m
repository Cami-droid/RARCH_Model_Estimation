% shutdown.m
% closing file to save outputs and result structs

% Definir el nombre de la carpeta donde se guardarán los resultados
resultsFolder = 'results';

% Crear la carpeta si no existe
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

% Definir el nombre del archivo para guardar las estructuras
outputFileName = fullfile(resultsFolder, 'project_results.mat');

% Guardar las estructuras outputs y results en el archivo
save(outputFileName, 'outputs', 'results');

% Mensaje indicando que el guardado fue exitoso
fprintf('Las estructuras outputs y results se han guardado exitosamente en %s\n', outputFileName);
