# RARCH_Model_Estimation
Proyecto m para generar parámetros y valores de verosimilitud para modelos multivariados rotados

incluye la estimación de 4 modelos multivariados en los cuales se calculas los parámetros ARCH(1) y GARCH(1),
las covarianzas condicionales sujetas a un conjunto de información disponible y las log-verosimilitudes.

Los 4 modelos incluidos en la estimación comparativa son:
-RBEKK
-OGARCH
-GOGARCH
-RDCC

Para correr el aplicativo debe ejecutarse deben cargarse a la consola matlab el archivo data_import_Table%n. script run_rarchs(n);n=2,4,5.

Para hacerlo correr con datos propios deben cargarse los datos al entorno matlab, generar un array(T,d) T es el tiempo y d el número de series. 
Las series deben ser de media 0. Luego ejecutar el archivo main.m

*************************************************************************************************************************************************

Project for generating parameters and log-likelihood values for rotated multivariate models.

It includes estimating 4 multivariate models in which the ARCH(1) and GARCH(!) parameters are computed, and the conditional covariance is subject 
to the available information set at this moment, and the log-likelihood values.

The 4 models included in the comparative estimation are:
-RBEKK
-OGARCH
-GOGARCH
-RDCC

To run the application, the file data_import_Table%n  must be loaded into the Matlab environment and then run the script run_rarchs(n); n=2,4,5.

If you wish to run the program with your data, you must load the data to the Matlab environment, and generate a (T,d) array where T is time and d 
the number of assets. The series should be demeaned. Then run the file main.m
