# modeloparametricogeodesia

Descripción de coords_EGM2008_2.txt
Este archivo es un conjunto de datos espaciales que representa el modelo del campo gravitatorio terrestre EGM2008, generado por el instituto gfz-potsdam. El archivo se estructura en dos partes principales: un bloque de encabezado con metadatos técnicos (que incluye el sistema de referencia GRS80, límites de latitud/longitud, tamaño de la cuadrícula, y valores estadísticos como medias y extremos) y un bloque de datos que comienza inmediatamente después de la línea end_of_head.   Descripción de las columnas
El bloque de datos del archivo está organizado en formato de cuadrícula (long_lat_value) y contiene las siguientes tres columnas:   longitude [deg.]: Indica la coordenada de longitud geográfica del punto evaluado, expresada en grados.   latitude [deg.]: Indica la coordenada de latitud geográfica del punto evaluado, expresada en grados.   geoid [meter]: Representa el valor funcional medido en la coordenada correspondiente, que en este caso es la altura u ondulación del geoide, expresada en metros.

El archivo "QGeoidCOL2023_3.txt" contiene una matriz de datos espaciales estructurada en tres columnas numéricas. De acuerdo con el formato estándar de los modelos geoidales, esta cuadrícula representa lo siguiente:   Latitud: La primera columna muestra coordenadas que descienden secuencialmente en intervalos regulares desde 12.800 hasta 11.633.   Longitud: La segunda columna indica las coordenadas longitudinales que avanzan de oeste a este, abarcando desde -79.299 hasta -66.716.   Ondulación Geoidal (N): La tercera columna refleja valores que varían entre 0.167 y -39.948, los cuales representan la altura del geoide (o cuasigeoide) respecto al elipsoide de referencia para cada par de coordenadas. 

"coords_QGEOID2023_3.txt". El conjunto de datos contiene 2318 registros estructurados en seis columnas separadas por punto y coma (;), correspondientes a una red de control geodésico y nivelación en Colombia.Estructura de los datos:FID: Identificador numérico secuencial del punto (0 a 2317).Altura_m_s: Altura física (ortométrica o normal) sobre el nivel medio del mar, en metros.Bloque_Aju: Agrupación o bloque de ajuste de la red (Circuito 1 y 2, Circuito 3 y 4, Circuito 5 y 6, Tramo 1).Latitud y Longitud: Coordenadas espaciales que trazan un recorrido longitudinal desde el suroccidente (aprox. 3.4°, -77.0°) hacia la región norte/Caribe (aprox. 11.5°, -72.4°).Alt_Elipso: Altura referida al elipsoide geodésico ($h$), en metros.


## Script `modelo_parametrico_zero_tide.py`

Solo usa la librería estándar de Python (3.8 o superior). Hace lo siguiente:

1. Lee los puntos. Pregunta qué columna es la latitud, la longitud, h (altura elipsoidal) y H (altura sobre el nivel del mar), y en qué sistema de mareas viene cada altura y cada malla.
2. Lleva todo a zero tide con las Ec. (1)-(3) de Guimarães et al. (2025). Lo que ya está en zero tide no se toca.
3. Depura los puntos con errores y los reporta en `puntos_rechazados.csv`:
   - R1: dato vacío o no numérico.
   - R2: coordenada fuera de rango.
   - R3: identificador o coordenadas duplicadas.
   - R4: h nula.
   - R5: \|h−H\| > 150 m.
   - R6: error grosero en h−H−N (por defecto, más de 5 σ robustas).
   - R7: punto fuera de la cobertura de la malla.
4. Exporta los puntos en zero tide a `coordenadas_zero_tide.csv`.
5. Interpola N en los puntos con los 6 métodos de `appQGeoidCOL2023_0.m`.
6. Ajusta el modelo de la Ec. (18), l = h−H−N = δH + x1(φ−φ0) + x2(λ−λ0)cosφ, con P = I, y calcula el R² ajustado de cada método. El mejor método es el de mayor R²adj.

```bash
python modelo_parametrico_zero_tide.py            # interactivo: pregunta columnas y sistemas de mareas
python modelo_parametrico_zero_tide.py --no-interactivo --marea-h zt --marea-H zt \
       --marea-malla zt --marea-malla zt          # sin preguntas (zt, mt o tf)
python modelo_parametrico_zero_tide.py --help     # resto de opciones
```

Si una malla no está en zero tide, se exporta convertida a `resultados/<malla>_zero_tide.txt`, con el mismo formato y el mismo orden del original. Por ejemplo, `coords_EGM2008_zero_tide.txt`. Para compararla nodo a nodo con una malla descargada en zero tide:

```bash
python comparar_mallas.py resultados/coords_EGM2008_zero_tide.txt EGM2008_ICGEM_zero_tide.gdf --csv diferencias.csv
```

Por defecto usa `coords_QGEOID2023.txt` como puntos, y `QGeoidCOL2023.txt` y `coords_EGM2008.txt` como mallas de N. Los resultados quedan en `resultados/`.
