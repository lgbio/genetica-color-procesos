- Se filtró dos traits: BerryC, PCTuberflesh. Se cambió background a blanco

- Se ejecutó Selección Genómica (GS) variando el número de marcadores así: 5, 10, 25, 50, 75, y 100.
- Los marcadores se obtuvieron de un proceso anterior de GWAS con las tres herramientas sin incluir Shesis
- Esos marcadores resultantes se les dió un puntaje con nuestra función de puntaje GSScore
- De esos marcadores se seleccionaron por cada trait HCL 200 marcadores, es decir para cada uno de los 28 traits (Tallo.C, Tallo.H, Tallo.L, ....) se seleccionaron los top 200.
- Y de allí se realizo la GS variando de 5, 10, y así sucesivamente.
- Igual, para la GS se toma todas las muestras o registros (600>) solo se limita el número de marcadores.
- De estos conjuntos se hace el proceso de GS de tomar un conjunto de entrenamiento (Geno, Phenos) y un conjunto de prueba (Geno, Phenos)
- La ejecución se realizó con una librería paralela de R que construí para todo estos procesos de GS y que por debajo utiliza la librería BGLR.
- Los algoritmos evaluados para GS fueron: BA BB BC BL BRR GBLUP RKHS. Excluí EGBLUP porque es muy costoso y da resultados menores o iguales a los otros.
- La validación cruzada se hizo para tres tiempos y tres folds (ver archivo .yml)

Finalmente,
- Te envio la función que construye el gráfico desde la tabla resultante de validación cruzada. Los datos y ejecución completa está en el punto 14x del pipeline que te envié la otra vez y que esta en el github.
- Para cambiar los nombres de los traits, solo edita el archivo de entrada "out-GS-crossval-varying-nmarkers.csv" y renombra cada trait con su nombre final (en inglés) y vuelve a ejecutar la función, debería crear los gráficos con los nombres nuevos. Igualmente los títulos de las gráficas los puedes cambiar dentro de la función.

Cualquier cosa me avisas.
