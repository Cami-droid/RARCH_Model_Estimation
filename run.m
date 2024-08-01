clear;clc;

% Solicitar al usuario que ingrese el número de tabla deseado
table_number = input('Ingrese el número de tabla que desea (2, 4, o 5): ');

% Validar la entrada del usuario
if ~ismember(table_number, [2, 4, 5])
    error('Número de tabla inválido. Por favor, ingrese 2, 4, o 5.');
end

% Definir el nombre del archivo de importación de datos según la tabla seleccionada
data_import_file = sprintf('data_import_Table%d', table_number);

if exist(data_import_file, 'file') == 2
    eval(data_import_file);
else
    error('El archivo %s no existe.', data_import_file);
end


% Ejecutar el script principal
main;
