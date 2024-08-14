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
