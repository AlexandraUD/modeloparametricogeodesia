function salida = appQGeoidCOL2023(modo, rutaArchivo, textoPuntos)
%APPQGEOIDCOL2023 Interfaz grafica -- Ejercicio No. 3, Sistemas de
%Posicionamiento Geodesico (2026-3): lectura, validacion, interpolacion y
%representacion grafica del modelo cuasigeoidal QGeoidCOL2023, calculo de
%la altura normal H^N = h - zeta y comparacion con GeoCol2004 y EGM2008.
%
%   Version complementada de completo6.m (Ejercicio No. 2, VEMOS2022). Se
%   conserva la arquitectura de esa app (archivo unico autocontenido,
%   estado en fig.UserData, bitacora, proyeccion Origen-Nacional, Kriging
%   Ordinario con variograma ajustado por Minimos Cuadrados Ponderados),
%   pero el cargue de datos, la validacion, los metodos y las salidas se
%   adaptan a una MALLA REGULAR de anomalia de altura (zeta).
%
%   USO:
%     >> appQGeoidCOL2023              % interfaz grafica. Si
%                                      % QGeoidCOL2023.txt esta en la misma
%                                      % carpeta que este archivo, se carga
%                                      % automaticamente; si no, use el
%                                      % boton "Cargar".
%     >> R = appQGeoidCOL2023('lote')  % modo por lotes (sin interfaz):
%                                      % valida la malla, evalua los seis
%                                      % metodos en los 5 puntos de control,
%                                      % compara con GeoCol2004 y EGM2008
%                                      % (valores_control_P1_P5.csv), controla
%                                      % con la nivelacion y ejecuta la
%                                      % validacion cruzada.
%     >> R = appQGeoidCOL2023('lote', ruta, {'P1: 4.6097, -74.0817, 2600'})
%     >> S = appQGeoidCOL2023('sensibilidad') % RMSE de validacion cruzada de
%                                      % IDW (p, k) y kriging (n, ventana)
%                                      % con los mismos nodos de control
%
%   FORMATO DE LOS PUNTOS (cuadro de texto o modo por lotes), uno por linea
%   o separados por ";":     [Nombre:] lat, lon [, h]
%     - lat, lon en grados decimales (punto decimal); h (altura
%       elipsoidal, m) es opcional. Sin h solo se reporta zeta.
%     - Las lineas que empiezan por "%" o "#" son comentarios.
%
%   ARCHIVO DE ENTRADA: QGeoidCOL2023.TXT, sin encabezado, tres columnas
%   numericas separadas por espacios: latitud (grados), longitud (grados),
%   anomalia de altura zeta (m). El archivo viene ordenado por filas de
%   latitud decreciente (N->S) y longitud creciente (W->E); la malla se
%   reconstruye por indices, por lo que el orden no afecta el resultado.
%
%   VALIDACIONES DE LA MALLA (cargarYValidarMalla): lineas numericas de 3
%   columnas, valores finitos, coordenadas en rango geografico valido,
%   nodos duplicados (y duplicados con zeta conflictiva), regularidad del
%   paso en latitud y longitud, dimensiones, nodos faltantes (malla
%   rectangular completa), resolucion espacial y cobertura.
%
%   METODOS (todos devuelven explicitamente nodos, pesos y contribuciones;
%   zeta = sum(w_i * zeta_i) y se verifica sum(w_i) = 1):
%     1. Vecino mas cercano  -- 1 nodo, w = 1 (distancia metrica).
%     2. Lineal por triangulos -- cada celda se divide SIEMPRE por la
%        diagonal SW-NE (division fija y reproducible; una triangulacion de
%        Delaunay sobre una malla cuadrada es ambigua porque los 4 vertices
%        son cociclicos). Pesos = coordenadas baricentricas (3 nodos).
%     3. Bilineal -- 4 nodos de la celda, pesos (1-u)(1-v), u(1-v), uv,
%        (1-u)v con u, v coordenadas normalizadas dentro de la celda.
%     4. Bicubica -- convolucion cubica de Keys (1981), a = -0.5 (el mismo
%        nucleo que interp2(...,'cubic') en mallas uniformes), 4x4 = 16
%        nodos, pesos w = wy(i)*wx(j). Requiere un nodo extra a cada lado:
%        en la celda de borde se informa "no aplicable" y NO se sustituye.
%     5. IDW -- k nodos mas cercanos, w_i ~ 1/d_i^p (Shepard, 1968).
%     6. Kriging ordinario -- semivariograma experimental local (ventana de
%        nodos alrededor del punto) ajustado por MCP (peso = numero de
%        pares) entre modelos esferico, exponencial, gaussiano y circular;
%        sistema de KO con los n vecinos mas cercanos. Reporta sigma_KO.
%   Los metodos 1-4 operan sobre los indices de la malla (lat/lon); IDW y
%   kriging requieren distancias euclidianas y proyectan a Origen-Nacional
%   MAGNA-SIRGAS (EPSG:9377, km), igual que en completo6.m. En IDW el
%   factor de escala de la proyeccion se cancela al normalizar los pesos.
%
%   COBERTURA: ningun metodo extrapola. Si el punto cae fuera del
%   rectangulo de la malla se reporta "fuera de cobertura" sin valor.
%
%   REDONDEO: los calculos se hacen en doble precision. En la presentacion
%   zeta y H^N se muestran a 0.001 m (1 mm, la resolucion con la que viene
%   zeta en el archivo), pesos a 1e-6 y contribuciones a 1e-4 m. Los CSV
%   exportados conservan la precision completa.
%
%   REPRODUCIBILIDAD: no hay pasos aleatorios sin semilla. La validacion
%   cruzada y la autoverificacion usan RandStream con semilla fija, sin
%   alterar el generador global de MATLAB.
%
%   Requiere MATLAB R2022a o posterior (uifigure, boxchart, clim).
%   No requiere Mapping Toolbox ni Statistics Toolbox.
%
%   Trabajo grupal -- Sistemas de Posicionamiento Geodesico, Ejercicio
%   No. 3 (2026-3).

if nargin >= 1 && (ischar(modo) || isstring(modo)) && strcmpi(modo, 'lote')
    if nargin < 2 || isempty(rutaArchivo)
        rutaArchivo = fullfile(fileparts(mfilename('fullpath')), 'QGeoidCOL2023.txt');
    end
    if nargin < 3 || isempty(textoPuntos)
        textoPuntos = puntosControlPredeterminados();
    end
    salida = ejecutarLote(rutaArchivo, textoPuntos);
    return
end
if nargin >= 1 && (ischar(modo) || isstring(modo)) && strcmpi(modo, 'sensibilidad')
    if nargin < 2 || isempty(rutaArchivo)
        rutaArchivo = fullfile(fileparts(mfilename('fullpath')), 'QGeoidCOL2023.txt');
    end
    salida = ejecutarSensibilidad(rutaArchivo);
    return
end

fig = uifigure('Name', 'QGeoidCOL2023 -- Anomalia de altura y altura normal', ...
    'Position', [40 40 1460 860]);

% ---- Estado compartido de la app (se guarda en fig.UserData) ----
app = struct();
app.M = []; app.reporte = []; app.verif = []; app.ruta = '';
app.datosListos = false;
app.puntos = table(); app.resultados = {}; app.metodoResultados = '';
app.comparacion = []; app.cv = []; app.modelosOut = table();

gl = uigridlayout(fig, [1 2]);
gl.ColumnWidth = {430, '1x'};

% =========================== PANEL IZQUIERDO ===========================
panelIzq = uipanel(gl, 'Title', 'Controles');
panelIzq.Layout.Row = 1; panelIzq.Layout.Column = 1;

glIzq = uigridlayout(panelIzq, [12 1]);
glIzq.RowHeight = {36, 30, 170, 20, 26, 84, 34, 100, 30, 30, 30, '1x'};
glIzq.Padding = [10 10 10 10]; glIzq.RowSpacing = 6;

app.lblInstruccion = uilabel(glIzq, 'Text', ...
    'Cargue la malla QGeoidCOL2023 en formato ".txt" (latitud, longitud, zeta)', ...
    'FontWeight', 'bold', 'WordWrap', 'on');
app.lblInstruccion.Layout.Row = 1;

app.btnCargar = uibutton(glIzq, 'Text', 'Cargar QGeoidCOL2023.TXT', ...
    'ButtonPushedFcn', @(src, evt) cbCargarArchivo(fig, ''));
app.btnCargar.Layout.Row = 2;

app.txtInfo = uitextarea(glIzq, 'Editable', 'off', 'Value', infoBase());
app.txtInfo.Layout.Row = 3;

app.lblMetodo = uilabel(glIzq, 'Text', 'Método de interpolación para "Calcular" y la vista local:');
app.lblMetodo.Layout.Row = 4;
app.ddMetodo = uidropdown(glIzq, 'Items', listaMetodos(), 'Value', 'Bilineal');
app.ddMetodo.Layout.Row = 5;

glPar = uigridlayout(glIzq, [2 4]);
glPar.Layout.Row = 6;
glPar.Padding = [0 0 0 0]; glPar.RowHeight = {40, 40}; glPar.RowSpacing = 4;
glPar.ColumnWidth = {'1.45x', '0.55x', '1.45x', '0.55x'};
l = uilabel(glPar, 'Text', 'IDW: Exponente p', 'WordWrap', 'on');
l.Tooltip = 'Potencia de la distancia: peso = 1/d^p. Mayor p = mas peso a los nodos cercanos.';
app.efP = uieditfield(glPar, 'numeric', 'Value', 2, 'Limits', [0.5 6]);
l = uilabel(glPar, 'Text', 'IDW: No. de nodos', 'WordWrap', 'on');
l.Tooltip = 'Cantidad de nodos mas cercanos usados por IDW.';
app.efK = uieditfield(glPar, 'numeric', 'Value', 16, 'Limits', [3 64], ...
    'RoundFractionalValues', 'on');
l = uilabel(glPar, 'Text', 'Kriging Ordinario: No. de vecinos', 'WordWrap', 'on');
l.Tooltip = 'Cantidad de nodos vecinos usados en cada prediccion de kriging.';
app.efNk = uieditfield(glPar, 'numeric', 'Value', 16, 'Limits', [4 64], ...
    'RoundFractionalValues', 'on');
l = uilabel(glPar, 'Text', 'Kriging Ordinario: tamaño de ventana (nodos)', 'WordWrap', 'on');
l.Tooltip = 'Semiancho de la ventana de busqueda alrededor del punto, en nodos de la malla.';
app.efVent = uieditfield(glPar, 'numeric', 'Value', 6, 'Limits', [3 12], ...
    'RoundFractionalValues', 'on');

app.lblPuntos = uilabel(glIzq, 'Text', ...
    'Punto(s): "[Nombre:] lat, lon [, h]" -- uno por linea o separados por ";" (h = altura elipsoidal, m)', ...
    'WordWrap', 'on');
app.lblPuntos.Layout.Row = 7;
app.txtPuntos = uitextarea(glIzq, 'Value', {'Ejemplo: 4.6097, -74.0817, 2600.000'});
app.txtPuntos.Layout.Row = 8;
app.txtPuntos.FontName = 'Consolas'; app.txtPuntos.FontSize = 9; % [NUEVO] columnas alineadas

glBot1 = uigridlayout(glIzq, [1 2]); glBot1.Layout.Row = 9; glBot1.Padding = [0 0 0 0];
app.btnCalcular = uibutton(glBot1, 'Text', 'Calcular (metodo activo)', 'Enable', 'off', ...
    'ButtonPushedFcn', @(src, evt) cbCalcular(fig));
app.btnComparar = uibutton(glBot1, 'Text', 'Comparar los 6 metodos', 'Enable', 'off', ...
    'ButtonPushedFcn', @(src, evt) cbCompararMetodos(fig));

glBot2 = uigridlayout(glIzq, [1 2]); glBot2.Layout.Row = 10; glBot2.Padding = [0 0 0 0];
app.btnPresets = uibutton(glBot2, 'Text', 'Cargar 5 puntos de control', ...
    'ButtonPushedFcn', @(src, evt) cbPuntosControl(fig));
app.btnValidar = uibutton(glBot2, 'Text', 'Validacion cruzada (malla)', 'Enable', 'off', ...
    'ButtonPushedFcn', @(src, evt) cbValidacionCruzada(fig));

app.btnExportar = uibutton(glIzq, 'Text', 'Exportar resultados (CSV, TXT y figuras PNG)', ...
    'Enable', 'off', 'ButtonPushedFcn', @(src, evt) cbExportar(fig, ''));
app.btnExportar.Layout.Row = 11;

app.txtLog = uitextarea(glIzq, 'Editable', 'off', ...
    'Value', {'BITACORA DE EJECUCION'; repmat('=', 1, 32); 'Esperando la carga de datos...'; ''}); % [MODIFICADO] encabezado
app.txtLog.Layout.Row = 12;
app.txtLog.FontName = 'Consolas'; app.txtLog.FontSize = 11; % [NUEVO] fuente monoespaciada

% =========================== PANEL DERECHO ============================
panelDer = uipanel(gl, 'Title', 'Resultados');
panelDer.Layout.Row = 1; panelDer.Layout.Column = 2;
glDer = uigridlayout(panelDer, [1 1]);
glDer.Padding = [6 6 6 6];

app.tabs = uitabgroup(glDer);

% --- Pestana 1: malla general y control de calidad ---
app.tabMalla = uitab(app.tabs, 'Title', 'Malla QGeoidCOL2023');
glMalla = uigridlayout(app.tabMalla, [1 2]);
glMalla.ColumnWidth = {'1.8x', '1x'};
app.axMalla = uiaxes(glMalla);
title(app.axMalla, 'Cargue los datos para ver la malla');
app.axHist = uiaxes(glMalla);
title(app.axHist, 'Histograma de \zeta');

% --- Pestana: control de calidad de la malla y autoverificacion ---
app.tabCalidad = uitab(app.tabs, 'Title', 'Control de calidad');
glCal = uigridlayout(app.tabCalidad, [1 1]);
app.tblCalidad = uitable(glCal);
app.tblCalidad.ColumnName = {'Grupo', 'Control', 'Resultado', 'Estado', 'Criterio / significado'};
app.tblCalidad.ColumnWidth = {120, 230, 230, 80, 'auto'};
app.tblCalidad.Data = cell(0, 5);
% [NUEVO] La pestana de control de calidad permanece oculta hasta cargar el .txt
app.tabCalidad.Parent = [];
% [FIN NUEVO]

% --- Pestana 2: consulta, nodos/pesos y vista local ---
app.tabConsulta = uitab(app.tabs, 'Title', 'Consulta y vista local');
glCons = uigridlayout(app.tabConsulta, [2 2]);
glCons.RowHeight = {'1.35x', '1x'}; glCons.ColumnWidth = {'1.25x', '1x'};
app.axLocal = uiaxes(glCons);
app.axLocal.Layout.Row = 1; app.axLocal.Layout.Column = 1;
title(app.axLocal, 'Vista local del campo interpolado');

glRes = uigridlayout(glCons, [2 2]);
glRes.Layout.Row = 1; glRes.Layout.Column = 2;
glRes.RowHeight = {26, '1x'}; glRes.ColumnWidth = {110, '1x'}; glRes.Padding = [0 0 0 0];
uilabel(glRes, 'Text', 'Punto en detalle:');
app.ddDetalle = uidropdown(glRes, 'Items', {'(sin calculos)'}, 'Enable', 'off', ...
    'ValueChangedFcn', @(src, evt) mostrarDetalle(fig, src.Value));
app.txtResumen = uitextarea(glRes, 'Editable', 'off', 'FontName', 'Consolas', ...
    'Value', {'Sin resultados.'});
app.txtResumen.Layout.Row = 2; app.txtResumen.Layout.Column = [1 2];

app.tblNodos = uitable(glCons);
app.tblNodos.Layout.Row = 2; app.tblNodos.Layout.Column = 1;
app.tblNodos.ColumnName = {'Nodo', 'Lat (grados)', 'Lon (grados)', 'zeta_i (m)', 'Peso w_i', ...
    'Contrib. w_i*zeta_i (m)', 'Distancia (km)'};
app.tblPuntos = uitable(glCons);
app.tblPuntos.Layout.Row = 2; app.tblPuntos.Layout.Column = 2;

% [NUEVO] --- Pestana: semivariograma / covariograma usado por el Kriging Ordinario ---
app.tabVario = uitab(app.tabs, 'Title', 'Variograma y covariograma (KO)');
glVar = uigridlayout(app.tabVario, [2 2]);
glVar.RowHeight = {'1.3x', '1x'};
app.axVario = uiaxes(glVar);
app.axVario.Layout.Row = 1; app.axVario.Layout.Column = 1;
title(app.axVario, 'Semivariograma ajustado');
app.axCovario = uiaxes(glVar);
app.axCovario.Layout.Row = 1; app.axCovario.Layout.Column = 2;
title(app.axCovario, 'Covariograma derivado');
app.tblVario = uitable(glVar);
app.tblVario.Layout.Row = 2; app.tblVario.Layout.Column = [1 2];
app.tblVario.ColumnName = {'Parametro', 'Valor', 'Unidad'};
app.tblVario.ColumnWidth = {330, 'auto', 90};
app.tblVario.Data = {'Estado', ['Seleccione "Kriging ordinario" como metodo y presione "Calcular" ' ...
    'para ver el variograma del punto en detalle.'], ''};
% [FIN NUEVO]

% --- Pestana 3: comparacion de metodos ---
app.tabMetodos = uitab(app.tabs, 'Title', 'Comparacion entre métodos de interpolacion');
glMet = uigridlayout(app.tabMetodos, [3 1]);
glMet.RowHeight = {28, '1x', '1x'};
glRef = uigridlayout(glMet, [1 3]); glRef.Padding = [0 0 0 0];
glRef.ColumnWidth = {170, 240, '1x'};
uilabel(glRef, 'Text', 'Metodo de referencia:');
app.ddReferencia = uidropdown(glRef, 'Items', listaMetodos(), 'Value', 'Bilineal', ...
    'ValueChangedFcn', @(src, evt) cbCompararMetodos(fig));
uilabel(glRef, 'Text', 'Dif. = zeta(método) - zeta(referencia), en m, cm y mm. La referencia tiene diferencia 0');
app.tblMetodos = uitable(glMet);
app.axMetodos = uiaxes(glMet);
title(app.axMetodos, 'Diferencias entre metodos');

% --- Pestana 4: validacion cruzada por decimacion ---
app.tabCV = uitab(app.tabs, 'Title', 'Validacion cruzada');
glCV = uigridlayout(app.tabCV, [3 1]);
glCV.RowHeight = {44, '0.8x', '1.3x'};
uilabel(glCV, 'WordWrap', 'on', 'Text', ['Se decima la malla (paso 10'', nodos de indice impar retirados), ' ...
    'se predicen los nodos retirados con cada metodo y se compara con su zeta real. ' ...
    'Muestra aleatoria con semilla fija (reproducible). Los errores son una cota ' ...
    'conservadora: el paso de prueba es el doble del paso real (5'').']);
app.tblCV = uitable(glCV);
glCVax = uigridlayout(glCV, [1 2]); glCVax.Padding = [0 0 0 0];
app.axCV = uiaxes(glCVax);
title(app.axCV, 'Distribucion de errores por metodo');
app.axCVbar = uiaxes(glCVax);
title(app.axCVbar, 'Indicadores de error por metodo');

% --- Pestana 5: comparacion entre modelos ---
app.tabModelos = uitab(app.tabs, 'Title', 'Comparación respecto GeoCol2004 y EGM2008');
glMod = uigridlayout(app.tabModelos, [4 1]);
glMod.RowHeight = {30, '0.8x', '0.8x', '1.25x'};
glModCtl = uigridlayout(glMod, [1 5]); glModCtl.Padding = [0 0 0 0];
glModCtl.ColumnWidth = {190, 180, 110, 250, '1x'};
app.btnLlenarMod = uibutton(glModCtl, 'Text', '1. Llenar con puntos actuales', ...
    'Enable', 'off', 'ButtonPushedFcn', @(src, evt) cbLlenarModelos(fig));
app.btnImportarMod = uibutton(glModCtl, 'Text', '2. Importar valores (CSV)', ...
    'Enable', 'off', 'ButtonPushedFcn', @(src, evt) cbImportarModelos(fig));
uilabel(glModCtl, 'Text', 'EGM2008 entrega:', 'HorizontalAlignment', 'right');
app.ddEGM = uidropdown(glModCtl, 'Items', ...
    {'zeta: anomalia de altura (cuasigeoide)', 'N: ondulacion geoidal (geoide)'}, ...
    'Value', 'N: ondulacion geoidal (geoide)');
app.btnCalcMod = uibutton(glModCtl, 'Text', '3. Calcular diferencias', ...
    'Enable', 'off', 'ButtonPushedFcn', @(src, evt) cbCalcularModelos(fig));
app.tblModIn = uitable(glMod, 'ColumnEditable', [false false false true false true true true]);
app.tblModIn.ColumnName = {'Punto', 'Lat', 'Lon', 'h (m)', 'zeta QGeoidCOL2023 (m)', ...
    'N GeoCol2004 (m) [editable]', 'EGM2008 (m) [editable]', 'H nivelada (m) [editable]'};
app.tblModOut = uitable(glMod);
glModAx = uigridlayout(glMod, [1 3]); glModAx.Padding = [0 0 0 0];
app.axDifGC = uiaxes(glModAx);  title(app.axDifGC, 'QGeoidCOL2023 - GeoCol2004');
app.axDifEGM = uiaxes(glModAx); title(app.axDifEGM, 'QGeoidCOL2023 - EGM2008');
app.axDifBar = uiaxes(glModAx); title(app.axDifBar, 'Diferencias por punto');

fig.UserData = app;

% Carga automatica si el archivo esta junto a este .m (no hay rutas
% absolutas propias del equipo; si no existe, se usa el boton "Cargar").
rutaDef = fullfile(fileparts(mfilename('fullpath')), 'QGeoidCOL2023.txt');
if isfile(rutaDef)
    cbCargarArchivo(fig, rutaDef);
end

if nargout > 0
    salida = fig;
end
end


% ============================================================
% CALLBACKS
% ============================================================

function cbCargarArchivo(fig, ruta)
app = fig.UserData;

if nargin < 2 || isempty(ruta)
    [archivo, carpeta] = uigetfile({'*.TXT;*.txt', 'Malla QGeoidCOL2023 (*.TXT)'}, ...
        'Seleccione el archivo QGeoidCOL2023.TXT');
    if isequal(archivo, 0)
        return
    end
    ruta = fullfile(carpeta, archivo);
end
[~, nom, ext] = fileparts(ruta);

app = agregarLog(app, sprintf('Importando "%s%s"...', nom, ext));
dlg = uiprogressdlg(fig, 'Title', 'Cargando malla', ...
    'Message', 'Leyendo, validando y autoverificando la implementacion...', 'Indeterminate', 'on');
try
    [M, rep] = cargarYValidarMalla(ruta);
    verif = autoverificarImplementacion(M);
catch ME
    close(dlg);
    uialert(fig, ME.message, 'Error al importar', 'Icon', 'error');
    app = agregarLog(app, ['ERROR al importar: ' ME.message]);
    fig.UserData = app;
    return
end
close(dlg);

app.M = M; app.reporte = rep; app.verif = verif; app.ruta = ruta;
app.datosListos = true;
app.puntos = table(); app.resultados = {}; app.comparacion = []; app.cv = [];
app.modelosOut = table();

set([app.btnCalcular, app.btnComparar, app.btnValidar, app.btnExportar, ...
    app.btnLlenarMod, app.btnImportarMod, app.btnCalcMod], 'Enable', 'on');

% [NUEVO] Se habilita (muestra) la pestana de control de calidad al cargar el .txt,
% en la segunda posicion (despues de la malla).
if isempty(app.tabCalidad.Parent)
    app.tabCalidad.Parent = app.tabs;
    ch = app.tabs.Children(:);
    app.tabs.Children = [ch(1); ch(end); ch(2:end-1)];
end
% [FIN NUEVO]

graficarMallaGeneral(app.axMalla, M, table());
graficarHistogramaZeta(app.axHist, M);

app.txtInfo.Value = [infoBase(); {'-----------------------------------------------'}; ...
    {sprintf('Malla: %d x %d nodos, paso %.4f'' x %.4f''', numel(M.lat), numel(M.lon), ...
        M.dLat*60, M.dLon*60)}; ...
    {sprintf('Cobertura: lat [%.6f, %.6f], lon [%.6f, %.6f]', ...
        M.lat(1), M.lat(end), M.lon(1), M.lon(end))}];

Tcal = tablaControlCalidad(rep, verif);
app.tblCalidad.Data = Tcal;
try
    removeStyle(app.tblCalidad);
catch
end
iRev = find(strcmp(Tcal(:, 4), 'REVISAR'));
if ~isempty(iRev)
    addStyle(app.tblCalidad, uistyle('BackgroundColor', [1 0.8 0.8]), 'row', iRev);
end
app = agregarLog(app, { ...
    sprintf('Malla cargada: %d x %d nodos, paso %.1f min.', rep.nLat, rep.nLon, rep.dLat*60); ...
    sprintf('Controles de calidad: %d OK, %d por revisar.', ...
        sum(strcmp(Tcal(:, 4), 'OK')), numel(iRev)); ...
    'Detalle completo en la pestana "Control de calidad".'});
app.tabs.SelectedTab = app.tabMalla; % [MODIFICADO] antes: app.tabCalidad
fig.UserData = app;
end


function cbCalcular(fig)
app = fig.UserData;
if ~app.datosListos
    uialert(fig, 'Primero cargue el archivo QGeoidCOL2023.TXT.', 'Sin datos');
    return
end
puntos = parsearPuntos(app.txtPuntos.Value);
if isempty(puntos)
    uialert(fig, 'No se reconocio ningun punto valido. Formato: [Nombre:] lat, lon [, h]', ...
        'Entrada invalida', 'Icon', 'warning');
    return
end

metodo = app.ddMetodo.Value;
par = leerParametros(app);
n = height(puntos);
resultados = cell(n, 1);
for k = 1:n
    resultados{k} = interpolarPunto(app.M, puntos.Lat(k), puntos.Lon(k), metodo, par, true);
end

app.puntos = puntos; app.resultados = resultados; app.metodoResultados = metodo;
Tres = tablaResultadosPuntos(puntos, resultados);
app.tblPuntos.Data = redondearParaPresentacion(Tres);
app.tblPuntos.ColumnName = {'Punto', 'Lat', 'Lon', 'Estado', 'zeta (m)', 'h (m)', ...
    'HN = h - zeta (m)', 'N. nodos', 'Suma pesos', 'sigma KO (m)', 'Observacion'};

app.ddDetalle.Items = cellstr(puntos.Punto);
app.ddDetalle.ItemsData = num2cell(1:n);
app.ddDetalle.Value = 1;
app.ddDetalle.Enable = 'on';

nFuera = sum(Tres.Estado ~= "Calculado");
% [MODIFICADO] formato de la bitacora
linCalc = {sprintf('Calculado con %s: %d punto(s), %d sin valor.', metodo, n, nFuera)};
if ~strcmp(textoParametros(par, metodo), 'sin parametros')
    linCalc{end+1, 1} = ['Parametros: ' textoParametros(par, metodo)];
end
app = agregarLog(app, linCalc);
for k = 1:n
    if ~resultados{k}.aplicable
        app = agregarLog(app, sprintf('  %s: %s', puntos.Punto(k), resultados{k}.motivo));
    end
end
fig.UserData = app;

graficarMallaGeneral(app.axMalla, app.M, puntos);
mostrarDetalle(fig, 1);
app.tabs.SelectedTab = app.tabConsulta;

if n == 1
    r = resultados{1};
    if r.aplicable
        msg = sprintf(['Punto %s: lat=%.6f, lon=%.6f\nMetodo: %s\n\n' ...
            'zeta = %.3f m\nh = %s\nHN = h - zeta = %s\n\nSuma de pesos = %.12f'], ...
            puntos.Punto(1), r.lat, r.lon, metodo, r.zeta, textoNum(puntos.h(1), '%.3f m'), ...
            textoNum(puntos.h(1) - r.zeta, '%.3f m'), r.sumaPesos);
    else
        msg = sprintf('Punto %s: %s', puntos.Punto(1), r.motivo);
    end
    uialert(fig, msg, 'Resultado de la interpolacion', 'Icon', 'info');
end
end


function mostrarDetalle(fig, k)
app = fig.UserData;
if isempty(app.resultados) || k > numel(app.resultados)
    return
end
r = app.resultados{k};
p = app.puntos(k, :);
par = leerParametros(app);

if r.aplicable
    T = r.nodos(:, {'Nodo', 'Lat', 'Lon', 'Zeta', 'Peso', 'Contrib', 'DistKm'});
    T.Lat = round(T.Lat, 6); T.Lon = round(T.Lon, 6); T.Zeta = round(T.Zeta, 3);
    T.Peso = round(T.Peso, 6); T.Contrib = round(T.Contrib, 4); T.DistKm = round(T.DistKm, 3);
    app.tblNodos.Data = T;
else
    app.tblNodos.Data = table();
end
app.txtResumen.Value = lineasResumenPunto(p, r, par);

dlg = uiprogressdlg(fig, 'Title', 'Vista local', 'Message', ...
    'Evaluando el campo local con el metodo activo...', 'Indeterminate', 'on');
campo = campoLocal(app.M, r, par);
close(dlg);
graficarVistaLocal(app.axLocal, app.M, r, campo, p.Punto);
% [NUEVO] variograma y covariograma del punto en detalle
graficarVariograma(app.axVario, app.axCovario, app.tblVario, r, par, p.Punto);
% [FIN NUEVO]
end


function cbCompararMetodos(fig)
app = fig.UserData;
if ~app.datosListos
    return
end
puntos = parsearPuntos(app.txtPuntos.Value);
if isempty(puntos)
    uialert(fig, 'No hay puntos validos para comparar.', 'Entrada invalida', 'Icon', 'warning');
    return
end
par = leerParametros(app);
ref = app.ddReferencia.Value;

dlg = uiprogressdlg(fig, 'Title', 'Comparando metodos', 'Message', ...
    sprintf('Evaluando 6 metodos en %d punto(s)...', height(puntos)), 'Indeterminate', 'on');
C = compararMetodos(app.M, puntos, par, ref);
close(dlg);
app.comparacion = C;

app.tblMetodos.Data = redondearParaPresentacion(C.tabla);
app.tblMetodos.ColumnName = {'Punto', 'Método de interpolacion', 'Zeta: anomalia (m)', ...
    'HN = h - zeta (m)', 'Dif. zeta vs ref. (m)', 'Dif. zeta vs ref. (cm)', ...
    'Dif. zeta vs ref. (mm)', 'Suma de pesos (=1)', 'Estado'};
graficarComparacionMetodos(app.axMetodos, C);

% [MODIFICADO] formato de la bitacora
lin = {sprintf('COMPARACION DE METODOS (ref.: %s)', ref); ...
    sprintf('Parametros: %s', textoParametros(par, 'todos'))};
for k = 1:height(puntos)
    if all(isnan(C.Z(k, :)))
        lin{end+1, 1} = sprintf('> %s: sin valor (fuera de cobertura)', puntos.Punto(k)); %#ok<AGROW>
    else
        lin{end+1, 1} = sprintf('> %s', puntos.Punto(k)); %#ok<AGROW>
        lin{end+1, 1} = sprintf('    dispersion %.1f mm | max|dif| %.1f mm (%s)', ...
            C.rangoMm(k), C.maxAbsMm(k), C.metodoMaxAbs(k)); %#ok<AGROW>
    end
end
app = agregarLog(app, lin);
app.tabs.SelectedTab = app.tabMetodos;
fig.UserData = app;
end


function cbValidacionCruzada(fig)
app = fig.UserData;
if ~app.datosListos
    return
end
par = leerParametros(app);
nControl = 400; semilla = 2023;
% [MODIFICADO] formato de la bitacora
app = agregarLog(app, {sprintf('Validacion cruzada por decimacion: %d nodos de control, semilla %d...', ...
    nControl, semilla); sprintf('Parametros: %s', textoParametros(par, 'todos'))});
dlg = uiprogressdlg(fig, 'Title', 'Validacion cruzada', 'Message', 'Iniciando...', 'Value', 0);
cv = validacionDecimacion(app.M, par, nControl, semilla, ...
    @(f, msg) actualizarProgreso(dlg, f, msg));
close(dlg);
app.cv = cv;

Tv = cv.tabla;
Tv{:, 2:end} = round(Tv{:, 2:end}, 3);
app.tblCV.Data = Tv;
app.tblCV.ColumnName = {'Metodo', 'Evaluados', 'No aplicables', 'Sesgo (mm)', 'RMSE (mm)', ...
    'MAE (mm)', 'P95 abs(e) (mm)', 'Max abs(e) (mm)'};
graficarValidacionCruzada(app.axCV, app.axCVbar, cv);

% [MODIFICADO] formato de la bitacora (tabla alineada, valores en mm)
lin = {sprintf('VALIDACION CRUZADA (malla decimada a %.0f'', %d controles)', cv.pasoPruebaMin, cv.nControl); ...
    sprintf('%-22s %7s %7s %7s', 'Metodo (mm)', 'RMSE', 'MAE', 'Max')};
for m = 1:height(cv.tabla)
    lin{end+1, 1} = sprintf('%-22s %7.1f %7.1f %7.1f', cv.tabla.Metodo(m), ...
        cv.tabla.RMSE_mm(m), cv.tabla.MAE_mm(m), cv.tabla.MaxAbs_mm(m)); %#ok<AGROW>
    if cv.tabla.NoAplicables(m) > 0
        lin{end+1, 1} = sprintf('   (%d nodos no aplicables)', cv.tabla.NoAplicables(m)); %#ok<AGROW>
    end
end
app = agregarLog(app, lin);
app.tabs.SelectedTab = app.tabCV;
fig.UserData = app;
end


function cbPuntosControl(fig)
app = fig.UserData;
app.txtPuntos.Value = puntosControlPredeterminados();
app = agregarLog(app, ['Puntos de control P1-P5 cargados (Guajira, Cordoba, Norte de Santander, ' ...
    'Bogota y Arauca) con su altura elipsoidal h.']);
fig.UserData = app;
end


function cbLlenarModelos(fig)
app = fig.UserData;
puntos = parsearPuntos(app.txtPuntos.Value);
if isempty(puntos)
    uialert(fig, 'No hay puntos validos en el cuadro de texto.', 'Sin puntos', 'Icon', 'warning');
    return
end
par = leerParametros(app);
metodo = app.ddMetodo.Value;
n = height(puntos);
z = nan(n, 1);
for k = 1:n
    r = interpolarPunto(app.M, puntos.Lat(k), puntos.Lon(k), metodo, par, false);
    z(k) = r.zeta;
end

% Se conservan los valores externos ya digitados para puntos con el mismo
% nombre, de modo que volver a llenar la tabla no los borra.
Nant = nan(n, 1); Eant = nan(n, 1); Hant = nan(n, 1);
Tprev = app.tblModIn.Data;
if istable(Tprev) && height(Tprev) > 0
    for k = 1:n
        m = find(Tprev.Punto == puntos.Punto(k), 1);
        if ~isempty(m)
            Nant(k) = Tprev.N_GeoCol2004(m); Eant(k) = Tprev.EGM2008(m);
            Hant(k) = Tprev.H_nivelada(m);
        end
    end
end
app.tblModIn.Data = table(puntos.Punto, puntos.Lat, puntos.Lon, puntos.h, z, Nant, Eant, Hant, ...
    'VariableNames', {'Punto', 'Lat', 'Lon', 'h', 'ZetaQ', 'N_GeoCol2004', 'EGM2008', 'H_nivelada'});
app.tblModIn.ColumnName = {'Punto', 'Lat', 'Lon', 'h (m)', ...
    sprintf('zeta QGeoidCOL2023 (m) [%s]', metodo), 'N GeoCol2004 (m) [editable]', ...
    'EGM2008 (m) [editable]', 'H nivelada (m) [editable]'};
app = agregarLog(app, sprintf(['Tabla de modelos llenada con %d punto(s); zeta QGeoidCOL2023 por %s. ' ...
    'Digite o importe los valores de GeoCol2004 y EGM2008.'], n, metodo));
app.tabs.SelectedTab = app.tabModelos;
fig.UserData = app;
end


function cbImportarModelos(fig)
app = fig.UserData;
Tin = app.tblModIn.Data;
if ~istable(Tin) || height(Tin) == 0
    uialert(fig, 'Primero use "1. Llenar con puntos actuales".', 'Tabla vacia', 'Icon', 'warning');
    return
end
[archivo, carpeta] = uigetfile({'*.csv;*.txt', 'Valores externos (*.csv, *.txt)'}, ...
    'Archivo con columnas Punto, GeoCol2004, EGM2008 [, H nivelada]');
if isequal(archivo, 0)
    return
end
try
    [Tin, nAsignados] = importarValoresExternos(fullfile(carpeta, archivo), Tin);
catch ME
    uialert(fig, ME.message, 'Error al importar', 'Icon', 'error');
    return
end
app.tblModIn.Data = Tin;
app = agregarLog(app, sprintf('Valores externos importados de "%s": %d punto(s) asignados.', ...
    archivo, nAsignados));
fig.UserData = app;
end


function cbCalcularModelos(fig)
app = fig.UserData;
Tin = app.tblModIn.Data;
if ~istable(Tin) || height(Tin) == 0
    uialert(fig, 'Primero use "1. Llenar con puntos actuales".', 'Tabla vacia', 'Icon', 'warning');
    return
end
egmEsZeta = startsWith(app.ddEGM.Value, 'zeta');
[Tout, est] = compararModelos(Tin);
app.modelosOut = Tout;

app.tblModOut.Data = redondearParaPresentacion(Tout);
app.tblModOut.ColumnName = {'Punto', 'Lat', 'Lon', 'h (m)', 'zeta Q (m)', 'N GC2004 (m)', ...
    'EGM2008 (m)', 'Q-GC (m)', 'Q-GC (cm)', 'Q-GC (mm)', 'Q-EGM (m)', 'Q-EGM (cm)', ...
    'Q-EGM (mm)', 'HN = h-zeta (m)', 'H = h-N GC (m)', 'H = h-EGM (m)', 'H nivelada (m)', ...
    'HN - Hniv (cm)', 'H(GC) - Hniv (cm)', 'H(EGM) - Hniv (cm)'};
graficarDiferenciasModelos(app.axDifGC, app.axDifEGM, app.axDifBar, app.M, Tout);

app = agregarLog(app, [{'-----------------------------------------------'}; ...
    lineasCompatibilidad(egmEsZeta); lineasEstadisticasModelos(est)]);
fig.UserData = app;
end


function cbExportar(fig, carpeta)
app = fig.UserData;
if nargin < 2 || isempty(carpeta)
    carpeta = uigetdir(pwd, 'Seleccione la carpeta para exportar resultados');
    if isequal(carpeta, 0)
        return
    end
end
sello = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
pref = fullfile(carpeta, ['QGeoidCOL2023_' sello '_']);
escritos = {};

writelines([lineasReporte(app.reporte); {''}; lineasVerificacion(app.verif)], ...
    [pref 'validacion_malla.txt']);
escritos{end+1} = 'validacion_malla.txt';
writelines(app.txtLog.Value, [pref 'bitacora.txt']);
escritos{end+1} = 'bitacora.txt';

if ~isempty(app.resultados)
    writetable(tablaResultadosPuntos(app.puntos, app.resultados), [pref 'resultados_' ...
        nombreArchivo(app.metodoResultados) '.csv']);
    escritos{end+1} = 'resultados_<metodo>.csv';
    Tn = table();
    for k = 1:numel(app.resultados)
        r = app.resultados{k};
        if r.aplicable
            Tk = r.nodos;
            Tk = addvars(Tk, repmat(app.puntos.Punto(k), height(Tk), 1), ...
                repmat(string(r.metodo), height(Tk), 1), 'Before', 1, ...
                'NewVariableNames', {'Punto', 'Metodo'});
            Tn = [Tn; Tk]; %#ok<AGROW>
        end
    end
    if height(Tn) > 0
        writetable(Tn, [pref 'nodos_pesos_contribuciones.csv']);
        escritos{end+1} = 'nodos_pesos_contribuciones.csv';
    end
end
if ~isempty(app.comparacion)
    writetable(app.comparacion.tabla, [pref 'comparacion_metodos.csv']);
    escritos{end+1} = 'comparacion_metodos.csv';
end
if ~isempty(app.cv)
    writetable(app.cv.tabla, [pref 'validacion_cruzada.csv']);
    escritos{end+1} = 'validacion_cruzada.csv';
end
if height(app.modelosOut) > 0
    writetable(app.modelosOut, [pref 'comparacion_modelos.csv']);
    escritos{end+1} = 'comparacion_modelos.csv';
end

ejes = {app.axMalla, 'fig_malla_general'; app.axHist, 'fig_histograma_zeta'; ...
    app.axLocal, 'fig_vista_local'; app.axMetodos, 'fig_comparacion_metodos'; ...
    app.axCV, 'fig_validacion_cruzada'; app.axCVbar, 'fig_validacion_cruzada_indicadores'; ...
    app.axDifGC, 'fig_dif_QGeoid_GeoCol2004'; ...
    app.axDifEGM, 'fig_dif_QGeoid_EGM2008'; app.axDifBar, 'fig_dif_modelos_barras'};
for k = 1:size(ejes, 1)
    if ~isempty(ejes{k, 1}.Children)
        exportgraphics(ejes{k, 1}, [pref ejes{k, 2} '.png'], 'Resolution', 300);
        escritos{end+1} = [ejes{k, 2} '.png']; %#ok<AGROW>
    end
end
app = agregarLog(app, [{sprintf('Exportados %d archivo(s) en "%s" con prefijo QGeoidCOL2023_%s_:', ...
    numel(escritos), carpeta, sello)}; strcat({'  '}, escritos(:))]);
fig.UserData = app;
end


% ============================================================
% MODO POR LOTES (sin interfaz grafica)
% ============================================================

function R = ejecutarLote(ruta, textoPuntos)
%EJECUTARLOTE Ejecuta el flujo completo sin interfaz y lo imprime en la
%ventana de comandos: validacion de la malla, autoverificacion, los seis
%metodos (con nodos y pesos) en los puntos dados, comparacion de metodos y
%validacion cruzada. Devuelve una estructura con todos los resultados.
[M, rep] = cargarYValidarMalla(ruta);
lin = lineasReporte(rep);
fprintf('%s\n', lin{:});
verif = autoverificarImplementacion(M);
lin = lineasVerificacion(verif);
fprintf('\n');
fprintf('%s\n', lin{:});

puntos = parsearPuntos(textoPuntos);
par = parametrosPredeterminados();
metodos = listaMetodos();
fprintf('\nParametros: %s\n', textoParametros(par, 'todos'));

detalle = cell(height(puntos), numel(metodos));
for k = 1:height(puntos)
    fprintf('\n==== %s  (lat=%.6f, lon=%.6f, h=%s) ====\n', puntos.Punto(k), ...
        puntos.Lat(k), puntos.Lon(k), textoNum(puntos.h(k), '%.3f m'));
    for m = 1:numel(metodos)
        r = interpolarPunto(M, puntos.Lat(k), puntos.Lon(k), metodos{m}, par, true);
        detalle{k, m} = r;
        if r.aplicable
            fprintf('  %-24s zeta=%9.4f m  HN=%s  nodos=%2d  sum(w)-1=%+.1e  %s\n', ...
                metodos{m}, r.zeta, textoNum(puntos.h(k) - r.zeta, '%.3f m'), ...
                height(r.nodos), r.sumaPesos - 1, r.info);
        else
            fprintf('  %-24s NO APLICABLE: %s\n', metodos{m}, r.motivo);
        end
    end
end

C = compararMetodos(M, puntos, par, 'Bilineal');
fprintf('\nCOMPARACION DE METODOS (referencia: Bilineal)\n');
disp(redondearParaPresentacion(C.tabla));

cv = validacionDecimacion(M, par, 400, 2023, []);
fprintf('\nVALIDACION CRUZADA (malla decimada, %d controles)\n', cv.nControl);
disp(cv.tabla);

% --- Comparacion con GeoCol2004 y EGM2008 y control con la nivelacion ---
% Los valores externos se leen de valores_control_P1_P5.csv (junto a este .m):
% columnas Punto, H_nivelada_m, GeoCol2004_N_m, EGM2008_N_m.
Tmod = table(); est = struct();
rutaExt = fullfile(fileparts(mfilename('fullpath')), 'valores_control_P1_P5.csv');
if isfile(rutaExt)
    n = height(puntos);
    z = nan(n, 1);
    for k = 1:n
        z(k) = interpolarPunto(M, puntos.Lat(k), puntos.Lon(k), 'Bilineal', par, false).zeta;
    end
    Tin = table(puntos.Punto, puntos.Lat, puntos.Lon, puntos.h, z, nan(n, 1), nan(n, 1), nan(n, 1), ...
        'VariableNames', {'Punto', 'Lat', 'Lon', 'h', 'ZetaQ', 'N_GeoCol2004', 'EGM2008', 'H_nivelada'});
    Tin = importarValoresExternos(rutaExt, Tin);
    [Tmod, est] = compararModelos(Tin);
    lin = [lineasCompatibilidad(false); lineasEstadisticasModelos(est)];
    fprintf('\n');
    fprintf('%s\n', lin{:});
    disp(redondearParaPresentacion(Tmod));
else
    fprintf('\n(No se encontro %s: se omite la comparacion con GeoCol2004 y EGM2008.)\n', rutaExt);
end

R = struct('M', M, 'reporte', rep, 'verificacion', verif, 'puntos', puntos, ...
    'parametros', par, 'detalle', {detalle}, 'comparacion', C, 'validacionCruzada', cv, ...
    'modelos', Tmod, 'estadisticasModelos', est);
end


function S = ejecutarSensibilidad(ruta)
%EJECUTARSENSIBILIDAD Sensibilidad de los parametros de IDW (p, k) y de
%kriging (n vecinos, semiancho de ventana) medida con la validacion
%cruzada por decimacion, usando los mismos 400 nodos de control (semilla
%2023) para que todas las configuraciones sean comparables.
M = cargarYValidarMalla(ruta);
nControl = 400; semilla = 2023;

pVals = [1 2 3 4]; kVals = [4 8 16 32];
filas = {};
for p = pVals
    for k = kVals
        par = parametrosPredeterminados(); par.pIDW = p; par.kIDW = k;
        cv = validacionDecimacion(M, par, nControl, semilla, [], {'IDW'});
        filas(end+1, :) = {"IDW", sprintf('p=%g, k=%d', p, k), cv.tabla.RMSE_mm, ...
            cv.tabla.MAE_mm, cv.tabla.MaxAbs_mm}; %#ok<AGROW>
        fprintf('IDW p=%g k=%2d  RMSE=%7.2f mm\n', p, k, cv.tabla.RMSE_mm);
    end
end
nVals = [8 16 32]; wVals = [4 6 9];
for n = nVals
    for w = wVals
        par = parametrosPredeterminados(); par.nKrig = n; par.ventanaKrig = w;
        cv = validacionDecimacion(M, par, nControl, semilla, [], {'Kriging ordinario'});
        filas(end+1, :) = {"Kriging ordinario", sprintf('n=%d, semiancho=%d', n, w), ...
            cv.tabla.RMSE_mm, cv.tabla.MAE_mm, cv.tabla.MaxAbs_mm}; %#ok<AGROW>
        fprintf('KO  n=%2d w=%d  RMSE=%7.2f mm\n', n, w, cv.tabla.RMSE_mm);
    end
end
S = cell2table(filas, 'VariableNames', {'Metodo', 'Parametros', 'RMSE_mm', 'MAE_mm', 'MaxAbs_mm'});
disp(S);
end


% ============================================================
% LECTURA Y VALIDACION DE LA MALLA
% ============================================================

function [M, rep] = cargarYValidarMalla(ruta)
%CARGARYVALIDARMALLA Lee QGeoidCOL2023.TXT (lat, lon, zeta), aplica los
%controles de calidad y reconstruye la malla regular.
%
%   M: estructura con lat (nLat x 1, ascendente), lon (nLon x 1,
%   ascendente), Z (nLat x nLon, NaN donde falte un nodo), dLat, dLon.
%   rep: reporte de validacion (ver lineasReporte).
if ~isfile(ruta)
    error('cargarYValidarMalla:archivoInexistente', 'No se encontro el archivo "%s".', ruta);
end
[~, nom, ext] = fileparts(ruta);
rep = struct();
rep.archivo = [nom ext];

texto = fileread(ruta);
lineas = strtrim(regexp(texto, '\r\n|\n|\r', 'split'));
lineas = lineas(:);
esVacia = cellfun(@isempty, lineas);
numLinea = find(~esVacia);
lineas = lineas(~esVacia);
n = numel(lineas);
rep.nLineas = n;

% --- 1. Datos numericos: cada linea debe tener exactamente 3 numeros ---
datos = nan(n, 3);
esNumerica = false(n, 1);
for k = 1:n
    [v, cnt, msg] = sscanf(lineas{k}, '%f');
    if cnt == 3 && isempty(msg)
        datos(k, :) = v.';
        esNumerica(k) = true;
    end
end
rep.nNoNumericas = sum(~esNumerica);
rep.lineasNoNumericas = numLinea(find(~esNumerica, 5)).';
datos = datos(esNumerica, :);

% --- 2. Valores finitos ---
noFinito = any(~isfinite(datos), 2);
rep.nNoFinitos = sum(noFinito);
datos(noFinito, :) = [];

% --- 3. Coordenadas dentro del rango geografico valido ---
fueraRango = abs(datos(:, 1)) > 90 | abs(datos(:, 2)) > 180;
rep.nFueraRango = sum(fueraRango);
datos(fueraRango, :) = [];
if size(datos, 1) < 16
    error('cargarYValidarMalla:sinDatos', ...
        'Quedan %d registros validos en "%s": no es posible reconstruir una malla.', ...
        size(datos, 1), rep.archivo);
end
lat = datos(:, 1); lon = datos(:, 2); z = datos(:, 3);

% Orden del archivo (solo informativo)
if lat(1) > lat(end)
    rep.orden = 'filas de latitud decreciente (N->S)';
else
    rep.orden = 'filas de latitud creciente (S->N)';
end
if lon(2) > lon(1)
    rep.orden = [rep.orden ', longitud creciente (W->E) en cada fila'];
else
    rep.orden = [rep.orden ', longitud no creciente en cada fila'];
end

% --- 4. Regularidad: paso constante en latitud y en longitud ---
% uniquetol agrupa valores que difieren solo por representacion binaria
% (p.ej. 12.800000000000001) con tolerancia de 1e-9 grados (~0.1 mm) y
% conserva valores reales del archivo, de modo que las desviaciones
% reportadas no quedan contaminadas por un redondeo previo.
latU = uniquetol(lat, 1e-9, 'DataScale', 1);
lonU = uniquetol(lon, 1e-9, 'DataScale', 1);
pasosLat = diff(latU); pasosLon = diff(lonU);
rep.dLatMed = median(pasosLat); rep.dLonMed = median(pasosLon);
rep.desvPasoLat = max(abs(pasosLat - rep.dLatMed));
rep.desvPasoLon = max(abs(pasosLon - rep.dLonMed));
tolPaso = 1e-7;
rep.esRegular = rep.desvPasoLat < tolPaso && rep.desvPasoLon < tolPaso;
if ~rep.esRegular
    error('cargarYValidarMalla:irregular', ['La malla no es regular: desviacion maxima ' ...
        'del paso %.3g grados en latitud y %.3g grados en longitud.'], ...
        rep.desvPasoLat, rep.desvPasoLon);
end

% --- 5. Dimensiones y malla nominal ---
nLat = numel(latU); nLon = numel(lonU);
rep.nLat = nLat; rep.nLon = nLon;
rep.nEsperados = nLat * nLon;
latVec = linspace(latU(1), latU(end), nLat).';
lonVec = linspace(lonU(1), lonU(end), nLon).';
dLat = (latVec(end) - latVec(1)) / (nLat - 1);
dLon = (lonVec(end) - lonVec(1)) / (nLon - 1);
iLat = round((lat - latVec(1)) / dLat) + 1;
iLon = round((lon - lonVec(1)) / dLon) + 1;
rep.desvNominal = max([abs(lat - latVec(iLat)); abs(lon - lonVec(iLon))]);

% --- 6. Duplicados (misma posicion) y duplicados con zeta distinta ---
lin = sub2ind([nLat nLon], iLat, iLon);
[linU, iPrimero] = unique(lin, 'stable');
rep.nDuplicados = numel(lin) - numel(linU);
zMax = accumarray(lin, z, [nLat*nLon 1], @max, NaN);
zMin = accumarray(lin, z, [nLat*nLon 1], @min, NaN);
rep.nDuplicadosConflicto = sum(zMax - zMin > 1e-6);

% --- 7. Reconstruccion y nodos faltantes ---
Z = nan(nLat, nLon);
Z(lin(iPrimero)) = z(iPrimero);
rep.nPresentes = sum(~isnan(Z(:)));
rep.nFaltantes = sum(isnan(Z(:)));
rep.esRectangularCompleta = rep.nFaltantes == 0;

% --- 8. Resolucion y cobertura (elipsoide GRS80) ---
a = 6378137.0; f = 1/298.257222101; e2 = 2*f - f^2;
phiM = deg2rad(mean(latVec([1 end])));
radioMeridiano = a*(1 - e2) / (1 - e2*sin(phiM)^2)^1.5;
radioNormal = a / sqrt(1 - e2*sin(phiM)^2);
rep.resLatKm = deg2rad(dLat) * radioMeridiano / 1000;
rep.resLonKm = deg2rad(dLon) * radioNormal * cos(phiM) / 1000;
rep.latLim = latVec([1 end]).'; rep.lonLim = lonVec([1 end]).';
rep.extensionNSkm = deg2rad(latVec(end) - latVec(1)) * radioMeridiano / 1000;
rep.extensionEOkm = deg2rad(lonVec(end) - lonVec(1)) * radioNormal * cos(phiM) / 1000;
rep.dLat = dLat; rep.dLon = dLon;

% --- 9. Estadisticas y plausibilidad de zeta ---
[rep.zMin, iMin] = min(Z(:)); [rep.zMax, iMax] = max(Z(:));
[r1, c1] = ind2sub(size(Z), iMin); [r2, c2] = ind2sub(size(Z), iMax);
rep.posMin = [latVec(r1), lonVec(c1)]; rep.posMax = [latVec(r2), lonVec(c2)];
rep.zMedia = mean(Z(:), 'omitnan'); rep.zDesv = std(Z(:), 'omitnan');
rep.nZetaImplausible = sum(abs(Z(:)) > 150); % |zeta| > 150 m no es plausible en la Tierra

M = struct('lat', latVec, 'lon', lonVec, 'Z', Z, 'dLat', dLat, 'dLon', dLon, ...
    'nombre', rep.archivo);
end


function lin = lineasReporte(rep)
%LINEASREPORTE Reporte de validacion de la malla como cell array de texto.
siNo = {'NO', 'SI'};
lin = { ...
    sprintf('CONTROL DE CALIDAD -- %s', rep.archivo); ...
    sprintf('  Lineas con contenido:                %d', rep.nLineas); ...
    sprintf('  Lineas no numericas o sin 3 columnas: %d %s', rep.nNoNumericas, ...
        textoLineas(rep.lineasNoNumericas)); ...
    sprintf('  Valores no finitos (NaN/Inf):        %d', rep.nNoFinitos); ...
    sprintf('  Coordenadas fuera de rango valido:   %d', rep.nFueraRango); ...
    sprintf('  Nodos duplicados:                    %d (con zeta distinta: %d)', ...
        rep.nDuplicados, rep.nDuplicadosConflicto); ...
    sprintf('  Dimensiones (lat x lon):             %d x %d = %d nodos esperados', ...
        rep.nLat, rep.nLon, rep.nEsperados); ...
    sprintf('  Nodos presentes / faltantes:         %d / %d', rep.nPresentes, rep.nFaltantes); ...
    sprintf('  Malla rectangular completa:          %s', siNo{rep.esRectangularCompleta + 1}); ...
    sprintf('  Malla regular:                       %s', siNo{rep.esRegular + 1}); ...
    sprintf('    paso lat = %.10f grados (%.4f''), desv. max. = %.1e grados', ...
        rep.dLat, rep.dLat*60, rep.desvPasoLat); ...
    sprintf('    paso lon = %.10f grados (%.4f''), desv. max. = %.1e grados', ...
        rep.dLon, rep.dLon*60, rep.desvPasoLon); ...
    sprintf('    desv. max. coordenadas vs malla nominal = %.1e grados', rep.desvNominal); ...
    sprintf('  Resolucion (lat media):              %.3f km (N-S) x %.3f km (E-O)', ...
        rep.resLatKm, rep.resLonKm); ...
    sprintf('  Cobertura: lat [%.6f, %.6f], lon [%.6f, %.6f]', ...
        rep.latLim(1), rep.latLim(2), rep.lonLim(1), rep.lonLim(2)); ...
    sprintf('  Extension aprox.:                    %.0f km (N-S) x %.0f km (E-O)', ...
        rep.extensionNSkm, rep.extensionEOkm); ...
    sprintf('  Orden del archivo: %s', rep.orden); ...
    sprintf('  zeta: min = %.3f m (%.4f, %.4f); max = %.3f m (%.4f, %.4f)', ...
        rep.zMin, rep.posMin, rep.zMax, rep.posMax); ...
    sprintf('        media = %.3f m; desv. estandar = %.3f m; |zeta| > 150 m: %d', ...
        rep.zMedia, rep.zDesv, rep.nZetaImplausible)};
end


function verif = autoverificarImplementacion(M)
%AUTOVERIFICARIMPLEMENTACION Contrasta las implementaciones propias contra
%interp2 de MATLAB en puntos aleatorios interiores (semilla fija) y
%verifica la suma de pesos de los seis metodos.
%   - Vecino mas cercano vs interp2 'nearest'
%   - Bilineal vs interp2 'linear'
%   - Bicubica vs interp2 'cubic' (convolucion cubica, a = -0.5)
%   Lineal por triangulos, IDW y kriging no tienen equivalente directo en
%   interp2; para ellos se verifica la suma de pesos y que el valor quede
%   dentro del rango de los nodos usados (IDW) o sea finito.
s = RandStream('mt19937ar', 'Seed', 1);
n = 300;
lat = M.lat(2) + (M.lat(end-2) - M.lat(2)) * rand(s, n, 1);
lon = M.lon(2) + (M.lon(end-2) - M.lon(2)) * rand(s, n, 1);
par = parametrosPredeterminados();

pares = {'Vecino mas cercano', 'nearest'; 'Bilineal', 'linear'; 'Bicubica', 'cubic'};
verif = struct('n', n, 'metodo', {pares(:, 1)}, 'maxDif', nan(3, 1));
for m = 1:size(pares, 1)
    zRef = interp2(M.lon, M.lat, M.Z, lon, lat, pares{m, 2});
    zPropio = nan(n, 1);
    for k = 1:n
        zPropio(k) = interpolarPunto(M, lat(k), lon(k), pares{m, 1}, par, false).zeta;
    end
    verif.maxDif(m) = max(abs(zPropio - zRef));
end

% Suma de pesos de los seis metodos (30 puntos, incluye kriging)
metodos = listaMetodos();
nW = 30;
verif.nPesos = nW;
verif.maxErrSumaPesos = nan(numel(metodos), 1);
for m = 1:numel(metodos)
    e = nan(nW, 1);
    for k = 1:nW
        r = interpolarPunto(M, lat(k), lon(k), metodos{m}, par, false);
        e(k) = abs(r.sumaPesos - 1);
    end
    verif.maxErrSumaPesos(m) = max(e);
end
verif.metodosPesos = metodos(:);
end


function lin = lineasVerificacion(verif)
lin = {sprintf('AUTOVERIFICACION DE LA IMPLEMENTACION (%d puntos interiores, semilla fija)', verif.n)};
for m = 1:numel(verif.metodo)
    lin{end+1, 1} = sprintf('  %-20s vs interp2: max|dif| = %.2e m', ...
        verif.metodo{m}, verif.maxDif(m)); %#ok<AGROW>
end
lin{end+1, 1} = sprintf('  Suma de pesos (%d puntos), max|1 - suma(w)| por metodo:', verif.nPesos);
for m = 1:numel(verif.metodosPesos)
    lin{end+1, 1} = sprintf('    %-24s %.1e', verif.metodosPesos{m}, verif.maxErrSumaPesos(m)); %#ok<AGROW>
end
end


% ============================================================
% NUCLEO DE INTERPOLACION: SEIS METODOS CON NODOS, PESOS Y CONTRIBUCIONES
% ============================================================

function metodos = listaMetodos()
metodos = {'Vecino mas cercano', 'Lineal por triangulos', 'Bilineal', 'Bicubica', ...
    'IDW', 'Kriging ordinario'};
end


function par = parametrosPredeterminados()
par = struct('pIDW', 2, 'kIDW', 16, 'nKrig', 16, 'ventanaKrig', 6, 'nLags', 12, ...
    'modelosKrig', {{'sph', 'exp', 'gau', 'cir'}}, 'modeloKrigFijo', []);
end


function par = leerParametros(app)
par = parametrosPredeterminados();
par.pIDW = app.efP.Value;
par.kIDW = round(app.efK.Value);
par.nKrig = round(app.efNk.Value);
par.ventanaKrig = round(app.efVent.Value);
end


function txt = textoParametros(par, metodo)
switch metodo
    case 'IDW'
        txt = sprintf('p = %g, k = %d', par.pIDW, par.kIDW);
    case 'Kriging ordinario'
        txt = sprintf('n = %d vecinos, semiancho de ventana = %d nodos', par.nKrig, par.ventanaKrig);
    case 'todos'
        txt = sprintf('IDW p = %g, k = %d; kriging n = %d, semiancho = %d nodos', ...
            par.pIDW, par.kIDW, par.nKrig, par.ventanaKrig);
    otherwise
        txt = 'sin parametros';
end
end


function res = interpolarPunto(M, lat0, lon0, metodo, par, conDetalle)
%INTERPOLARPUNTO Interpola zeta en (lat0, lon0) con el metodo indicado.
%Cada metodo devuelve los indices de los nodos usados y sus pesos; zeta
%se calcula aqui como sum(w_i * zeta_i), igual para todos los metodos. Con
%conDetalle = true se arma la tabla de nodos (coordenadas, zeta, peso,
%contribucion y distancia) para el reporte.
%
%   Nunca se extrapola: si el punto esta fuera de la malla, o el metodo no
%   tiene soporte espacial suficiente, se devuelve aplicable = false y el
%   motivo, sin sustituir silenciosamente por otro metodo.
if nargin < 6
    conDetalle = true;
end
res = struct('metodo', metodo, 'lat', lat0, 'lon', lon0, 'zeta', NaN, ...
    'aplicable', false, 'motivo', '', 'nodos', table(), 'sumaPesos', NaN, ...
    'sigma', NaN, 'geom', crearGeom('', [], [], [], []), 'info', '', 'modeloKrig', []);

if ~dentroDeCobertura(M, lat0, lon0)
    res.motivo = sprintf(['Fuera de la cobertura de la malla (lat [%.4f, %.4f], ' ...
        'lon [%.4f, %.4f]): no se extrapola.'], M.lat(1), M.lat(end), M.lon(1), M.lon(end));
    return
end

extra = struct('sigma', NaN, 'info', '', 'modeloKrig', []);
switch metodo
    case 'Vecino mas cercano'
        [iLat, iLon, w, geom, motivo] = pesosVecinoMasCercano(M, lat0, lon0);
    case 'Lineal por triangulos'
        [iLat, iLon, w, geom, motivo] = pesosTriangulos(M, lat0, lon0);
    case 'Bilineal'
        [iLat, iLon, w, geom, motivo] = pesosBilineal(M, lat0, lon0);
    case 'Bicubica'
        [iLat, iLon, w, geom, motivo] = pesosBicubica(M, lat0, lon0);
    case 'IDW'
        [iLat, iLon, w, geom, motivo] = pesosIDW(M, lat0, lon0, par.pIDW, par.kIDW);
    case 'Kriging ordinario'
        [iLat, iLon, w, geom, motivo, extra] = pesosKriging(M, lat0, lon0, par);
    otherwise
        error('interpolarPunto:metodoDesconocido', 'Metodo no soportado: %s', metodo);
end
res.geom = geom;
res.sigma = extra.sigma; res.info = extra.info; res.modeloKrig = extra.modeloKrig;
if ~isempty(motivo)
    res.motivo = motivo;
    return
end

iLat = iLat(:); iLon = iLon(:); w = w(:);
z = M.Z(sub2ind(size(M.Z), iLat, iLon));
if any(isnan(z(w ~= 0)))
    res.motivo = 'Soporte incompleto: falta en la malla al menos un nodo con peso no nulo.';
    return
end
zc = z; zc(isnan(zc)) = 0; % nodos faltantes con peso exactamente cero
res.zeta = sum(w .* zc);
res.sumaPesos = sum(w);
res.aplicable = true;
if isempty(res.info)
    res.info = sprintf('%d nodo(s)', numel(w));
end
if conDetalle
    d = distanciaPlana(M.lat(iLat), M.lon(iLon), lat0, lon0) / 1000;
    res.nodos = table((1:numel(w)).', M.lat(iLat), M.lon(iLon), z, w, w .* z, d, iLat, iLon, ...
        'VariableNames', {'Nodo', 'Lat', 'Lon', 'Zeta', 'Peso', 'Contrib', 'DistKm', 'iLat', 'iLon'});
end
end


function g = crearGeom(tipo, lat, lon, celdaLat, celdaLon)
%CREARGEOM Geometria de soporte para dibujar en la vista local: 'vecino'
%(segmento punto-nodo), 'poligono' (triangulo, celda o ventana 4x4) o
%'radial' (segmentos del punto a cada nodo, IDW y kriging).
g = struct('tipo', tipo, 'lat', lat, 'lon', lon, 'celdaLat', celdaLat, 'celdaLon', celdaLon);
end


function [iLat, iLon, w, geom, motivo] = pesosVecinoMasCercano(M, lat0, lon0)
%PESOSVECINOMASCERCANO Nodo mas cercano (distancia metrica en Origen-
%Nacional). En una malla regular el nodo mas cercano siempre es una de las
%4 esquinas de la celda que contiene al punto.
motivo = '';
[i, j] = ubicarCelda(M, lat0, lon0);
cI = [i; i; i+1; i+1]; cJ = [j; j+1; j+1; j];
d = distanciaPlana(M.lat(cI), M.lon(cJ), lat0, lon0);
[~, k] = min(d); % en empate exacto se toma la primera esquina (SW, SE, NE, NW)
iLat = cI(k); iLon = cJ(k); w = 1;
geom = crearGeom('vecino', [lat0; M.lat(iLat)], [lon0; M.lon(iLon)], ...
    M.lat([i i i+1 i+1 i]), M.lon([j j+1 j+1 j j]));
end


function [iLat, iLon, w, geom, motivo] = pesosTriangulos(M, lat0, lon0)
%PESOSTRIANGULOS Interpolacion lineal en el triangulo que contiene al
%punto. La celda se divide por la diagonal SW-NE:
%   u >= v: triangulo SW-SE-NE, pesos [1-u, u-v, v]
%   u <  v: triangulo SW-NE-NW, pesos [1-v, u, v-u]
%(coordenadas baricentricas; en el triangulo la superficie es un plano).
motivo = '';
[i, j, u, v] = ubicarCelda(M, lat0, lon0);
if u >= v
    iLat = [i; i; i+1]; iLon = [j; j+1; j+1]; w = [1-u; u-v; v];
else
    iLat = [i; i+1; i+1]; iLon = [j; j+1; j]; w = [1-v; u; v-u];
end
geom = crearGeom('poligono', M.lat([iLat; iLat(1)]), M.lon([iLon; iLon(1)]), ...
    M.lat([i i i+1 i+1 i i+1]), M.lon([j j+1 j+1 j j j+1]));
end


function [iLat, iLon, w, geom, motivo] = pesosBilineal(M, lat0, lon0)
%PESOSBILINEAL Interpolacion bilineal en la celda (esquinas SW, SE, NE, NW).
motivo = '';
[i, j, u, v] = ubicarCelda(M, lat0, lon0);
iLat = [i; i; i+1; i+1]; iLon = [j; j+1; j+1; j];
w = [(1-u)*(1-v); u*(1-v); u*v; (1-u)*v];
geom = crearGeom('poligono', M.lat([iLat; iLat(1)]), M.lon([iLon; iLon(1)]), [], []);
end


function [iLat, iLon, w, geom, motivo] = pesosBicubica(M, lat0, lon0)
%PESOSBICUBICA Convolucion cubica de Keys (1981), a = -0.5, sobre los 4x4
%nodos alrededor de la celda. Pesos separables w(i,j) = wy(i)*wx(j); cada
%vector de pesos 1D suma 1, por lo que la suma total tambien es 1. Algunos
%pesos son negativos (lobulos del nucleo): es lo esperado.
motivo = ''; iLat = []; iLon = []; w = [];
[i, j, u, v] = ubicarCelda(M, lat0, lon0);
geomCelda = crearGeom('poligono', [], [], M.lat([i i i+1 i+1 i]), M.lon([j j+1 j+1 j j]));
if i - 1 < 1 || i + 2 > numel(M.lat) || j - 1 < 1 || j + 2 > numel(M.lon)
    geom = geomCelda;
    motivo = sprintf(['No aplicable: soporte 4x4 incompleto (celda de borde: fila %d de %d, ' ...
        'columna %d de %d). La convolucion cubica necesita un nodo adicional a cada lado; ' ...
        'no se sustituye por otro metodo.'], i, numel(M.lat) - 1, j, numel(M.lon) - 1);
    return
end
wx = nucleoCubico([u + 1, u, 1 - u, 2 - u]); % nodos lon j-1, j, j+1, j+2
wy = nucleoCubico([v + 1, v, 1 - v, 2 - v]); % nodos lat i-1, i, i+1, i+2
[JJ, II] = meshgrid(j-1:j+2, i-1:i+2);
W = wy(:) * wx(:).';
iLat = II(:); iLon = JJ(:); w = W(:);
geom = crearGeom('poligono', M.lat([i-1 i-1 i+2 i+2 i-1]), M.lon([j-1 j+2 j+2 j-1 j-1]), ...
    geomCelda.celdaLat, geomCelda.celdaLon);
end


function w = nucleoCubico(s)
%NUCLEOCUBICO Nucleo de convolucion cubica de Keys (1981) con a = -0.5.
a = -0.5;
s = abs(s);
w = zeros(size(s));
m1 = s <= 1;
m2 = s > 1 & s < 2;
w(m1) = (a + 2)*s(m1).^3 - (a + 3)*s(m1).^2 + 1;
w(m2) = a*s(m2).^3 - 5*a*s(m2).^2 + 8*a*s(m2) - 4*a;
end


function [iLat, iLon, w, geom, motivo] = pesosIDW(M, lat0, lon0, p, k)
%PESOSIDW Distancia Inversa Ponderada con los k nodos mas cercanos
%(distancias en Origen-Nacional). Si el punto coincide con un nodo (< 1 mm)
%ese nodo recibe peso 1. Cerca del borde de la malla los k nodos quedan
%todos hacia el interior (efecto de borde).
motivo = '';
[i, j] = ubicarCelda(M, lat0, lon0);
r = ceil(sqrt(k)/2) + 1;
[iLat, iLon] = ventanaIndices(M, i, j, r);
d = distanciaPlana(M.lat(iLat), M.lon(iLon), lat0, lon0);
[dOrd, orden] = sort(d);
kk = min(k, numel(d));
iLat = iLat(orden(1:kk)); iLon = iLon(orden(1:kk)); dOrd = dOrd(1:kk);
if dOrd(1) < 1e-3
    w = zeros(kk, 1); w(1) = 1;
else
    w = 1 ./ dOrd.^p;
    w = w / sum(w);
end
geom = crearGeom('radial', M.lat(iLat), M.lon(iLon), [], []);
end


function [iLat, iLon, w, geom, motivo, extra] = pesosKriging(M, lat0, lon0, par)
%PESOSKRIGING Kriging ordinario local. El semivariograma se ajusta por MCP
%con los nodos de una ventana de (2*semiancho)^2 nodos alrededor del punto
%(estacionariedad local: zeta tiene tendencias regionales de decenas de
%metros, por lo que un variograma global no representa la variabilidad a
%la escala de la interpolacion). Si par.modeloKrigFijo no esta vacio, se
%reutiliza ese modelo (vista local: un solo ajuste para todo el campo).
motivo = ''; iLat = []; iLon = []; w = [];
extra = struct('sigma', NaN, 'info', '', 'modeloKrig', []);
[i, j] = ubicarCelda(M, lat0, lon0);
semiancho = max(par.ventanaKrig, ceil(sqrt(par.nKrig)/2) + 1);
[vLat, vLon] = ventanaIndices(M, i, j, semiancho);
geom = crearGeom('radial', [], [], [], []);

[x, y] = proyectarOrigenNacional(M.lat(vLat), M.lon(vLon));
xk = x / 1000; yk = y / 1000; % km
z = M.Z(sub2ind(size(M.Z), vLat, vLon));
[x0, y0] = proyectarOrigenNacional(lat0, lon0);
x0 = x0 / 1000; y0 = y0 / 1000;

if isempty(par.modeloKrigFijo)
    try
        modelo = ajustarVariogramaMCP(xk, yk, z, par.nLags, par.modelosKrig);
    catch ME
        motivo = ['No aplicable: no fue posible ajustar el variograma local (' ME.message ').'];
        return
    end
else
    modelo = par.modeloKrigFijo;
end

[lambda, mu, varPred, rc, idx] = kriegeSistema(xk, yk, x0, y0, modelo, par.nKrig);
iLat = vLat(idx); iLon = vLon(idx); w = lambda;
geom = crearGeom('radial', M.lat(iLat), M.lon(iLon), [], []);
extra.modeloKrig = modelo;
extra.sigma = sqrt(max(varPred, 0));
extra.info = sprintf(['variograma %s (MCP, %d nodos): pepita=%.3g m2, meseta parcial=%.3g m2, ' ...
    'rango=%.1f km; mu=%.3g; sigma_KO=%.4f m; rcond=%.1e'], nombreModelo(modelo.tipo), ...
    modelo.nNodos, modelo.nugget, modelo.sill, modelo.rango, mu, extra.sigma, rc);
if rc < 1e-14
    extra.info = [extra.info ' (ADVERTENCIA: sistema mal condicionado)'];
end
end


function [i, j, u, v] = ubicarCelda(M, lat0, lon0)
%UBICARCELDA Celda [lat(i), lat(i+1)] x [lon(j), lon(j+1)] que contiene al
%punto y sus coordenadas normalizadas u (lon) y v (lat) en [0, 1].
i = ubicarCelda1D(M.lat, lat0);
j = ubicarCelda1D(M.lon, lon0);
v = min(max((lat0 - M.lat(i)) / (M.lat(i+1) - M.lat(i)), 0), 1);
u = min(max((lon0 - M.lon(j)) / (M.lon(j+1) - M.lon(j)), 0), 1);
end


function [idxLo, ok] = ubicarCelda1D(malla, val)
%UBICARCELDA1D Ubica, en un vector ordenado "malla", el indice idxLo tal
%que malla(idxLo) <= val <= malla(idxLo+1) (con idxLo entre 1 y
%numel(malla)-1). Si val coincide (con tolerancia) con el valor maximo de
%la malla, se usa la ultima celda (idxLo=numel(malla)-1), de forma que el
%peso de esa esquina resulte 1 y el de la esquina "siguiente" (que no
%existe) resulte 0. Devuelve ok=false si val cae fuera de
%[malla(1), malla(end)].
tol = 1e-9;
if numel(malla) < 2 || val < malla(1) - tol || val > malla(end) + tol
    idxLo = NaN; ok = false;
    return
end
idxLo = find(malla <= val + tol, 1, 'last');
if idxLo >= numel(malla)
    idxLo = numel(malla) - 1;
end
ok = true;
end


function [iLat, iLon] = ventanaIndices(M, i, j, r)
%VENTANAINDICES Indices de los nodos (no faltantes) de una ventana de
%2r x 2r nodos centrada en la celda (i, j), recortada al borde de la malla.
ii = max(1, i - r + 1):min(numel(M.lat), i + r);
jj = max(1, j - r + 1):min(numel(M.lon), j + r);
[JJ, II] = meshgrid(jj, ii);
iLat = II(:); iLon = JJ(:);
ok = ~isnan(M.Z(sub2ind(size(M.Z), iLat, iLon)));
iLat = iLat(ok); iLon = iLon(ok);
end


function dentro = dentroDeCobertura(M, lat0, lon0)
%DENTRODECOBERTURA Indica si (lat0, lon0) esta dentro del rectangulo de la
%malla (limites incluidos). Fuera de el ningun metodo se evalua.
tol = 1e-9;
dentro = lat0 >= M.lat(1) - tol & lat0 <= M.lat(end) + tol & ...
    lon0 >= M.lon(1) - tol & lon0 <= M.lon(end) + tol;
end


function d = distanciaPlana(lat, lon, lat0, lon0)
%DISTANCIAPLANA Distancia (m) en la proyeccion Origen-Nacional entre los
%nodos (lat, lon) y el punto (lat0, lon0).
[x, y] = proyectarOrigenNacional(lat(:), lon(:));
[x0, y0] = proyectarOrigenNacional(lat0, lon0);
d = hypot(x - x0, y - y0);
end


function [x, y] = proyectarOrigenNacional(lat, lon)
%PROYECTARORIGENNACIONAL Proyeccion MAGNA-SIRGAS / Origen-Nacional
%(EPSG:9377), formulas cerradas de Snyder (1987) sobre el elipsoide GRS80.
a = 6378137.0;
invF = 298.257222101;
f = 1 / invF;
e2 = 2*f - f^2;
ep2 = e2 / (1 - e2);

lat0 = deg2rad(4.0);
lon0 = deg2rad(-73.0);
k0 = 0.9992;
FE = 5000000.0;
FN = 2000000.0;

phi = deg2rad(lat(:));
lambda = deg2rad(lon(:));

M = @(p) a * ( ...
    (1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) .* p ...
    - (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) .* sin(2*p) ...
    + (15*e2^2/256 + 45*e2^3/1024) .* sin(4*p) ...
    - (35*e2^3/3072) .* sin(6*p));

M0 = M(lat0);
Mp = M(phi);

N = a ./ sqrt(1 - e2 .* sin(phi).^2);
Tt = tan(phi).^2;
C = ep2 .* cos(phi).^2;
A = (lambda - lon0) .* cos(phi);

x = FE + k0 .* N .* (A + (1 - Tt + C).*A.^3/6 ...
    + (5 - 18*Tt + Tt.^2 + 72*C - 58*ep2).*A.^5/120);

y = FN + k0 .* (Mp - M0 + N.*tan(phi).*(A.^2/2 ...
    + (5 - Tt + 9*C + 4*C.^2).*A.^4/24 ...
    + (61 - 58*Tt + Tt.^2 + 600*C - 330*ep2).*A.^6/720));

x = reshape(x, size(lat));
y = reshape(y, size(lat));
end


% ============================================================
% KRIGING ORDINARIO AJUSTADO POR MINIMOS CUADRADOS PONDERADOS (MCP)
% (Adaptado de completo6.m: ahora recibe coordenadas planas en km y la
% variable directamente, usa todos los pares de la ventana local -- sin
% muestreo aleatorio, para que el resultado sea reproducible -- y fija
% gamma(0) = 0 para que el kriging sea exacto en los nodos.)
% ============================================================

function [modeloElegido, tablaComparativa] = ajustarVariogramaMCP(x, y, z, nLags, tipos)
%AJUSTARVARIOGRAMAMCP Ajusta el semivariograma experimental de z por
%Minimos Cuadrados Ponderados (peso = numero de pares por rezago;
%equivalente a weights="npairs" en geoR::variofit o fit.method=2 en
%gstat::fit.variogram), probando varios modelos y quedandose con el de
%menor SSR ponderado. x, y en km.
x = x(:); y = y(:); z = z(:);
n = numel(z);
if n < 6
    error('ajustarVariogramaMCP:pocosDatos', 'solo %d nodos en la ventana', n);
end

D = hypot(x - x.', y - y.');
S = 0.5 * (z - z.').^2;
mascara = triu(true(n), 1);
dist = D(mascara); semiv = S(mascara);

lagMax = max(dist) / 2;
edges = linspace(0, lagMax, nLags + 1);
h = nan(nLags, 1); gammaExp = nan(nLags, 1); nPares = zeros(nLags, 1);
for k = 1:nLags
    enBin = dist > edges(k) & dist <= edges(k+1);
    nPares(k) = sum(enBin);
    if nPares(k) > 0
        h(k) = mean(dist(enBin));
        gammaExp(k) = mean(semiv(enBin));
    end
end
val = nPares > 0;
h = h(val); gammaExp = gammaExp(val); nPares = nPares(val);
if numel(h) < 3
    error('ajustarVariogramaMCP:pocosRezagos', 'menos de 3 rezagos con pares');
end

% El rango se reparametriza con una sigmoide acotada a (0, lagMax] (ver
% completo6.m: sin esta cota fminsearch puede divergir a rangos absurdos
% cuando el semivariograma experimental no se estabiliza en la ventana).
res = struct('tipo', {}, 'nugget', {}, 'sill', {}, 'rango', {}, 'SSRw', {});
ini = [sqrt(max(min(gammaExp)*0.3, eps)), sqrt(max(max(gammaExp) - min(gammaExp), eps)), 0];
opts = optimset('Display', 'off');
for i = 1:numel(tipos)
    obj = @(p) sum(nPares .* (gammaExp - modeloVariograma(h, p(1)^2, p(2)^2, ...
        lagMax/(1+exp(-p(3))) + eps, tipos{i})).^2);
    p = fminsearch(obj, ini, opts);
    res(i).tipo = tipos{i};
    res(i).nugget = p(1)^2;
    res(i).sill = p(2)^2;
    res(i).rango = lagMax/(1+exp(-p(3))) + eps;
    res(i).SSRw = obj(p);
end
tablaComparativa = struct2table(res);

[~, iMejor] = min(tablaComparativa.SSRw);
fila = tablaComparativa(iMejor, :);
modeloElegido = struct('tipo', fila.tipo{1}, 'nugget', fila.nugget, ...
    'sill', fila.sill, 'rango', fila.rango, 'SSRw', fila.SSRw, 'lagMax', lagMax, ...
    'nNodos', n, 'hExp', h, 'gammaExp', gammaExp, 'nPares', nPares);
end


function gamma = modeloVariograma(h, nugget, sill, rango, tipo)
%MODELOVARIOGRAMA Modelos teoricos de semivariograma candidatos. Por
%definicion gamma(0) = 0 (la pepita es el salto en el origen, h -> 0+).
h = h(:);
switch tipo
    case 'sph'
        ratio = min(h/rango, 1);
        gamma = nugget + sill * (1.5*ratio - 0.5*ratio.^3);
    case 'exp'
        gamma = nugget + sill * (1 - exp(-h/rango));
    case 'gau'
        gamma = nugget + sill * (1 - exp(-(h/rango).^2));
    case 'cir'
        ratio = min(h/rango, 1);
        gamma = nugget + sill * (1 - (2/pi)*(acos(ratio) - ratio.*sqrt(1 - ratio.^2)));
    otherwise
        error('modeloVariograma:tipoNoSoportado', 'Modelo no soportado: %s', tipo);
end
gamma(h <= 1e-12) = 0;
end


function txt = nombreModelo(tipo)
switch tipo
    case 'sph', txt = 'esferico';
    case 'exp', txt = 'exponencial';
    case 'gau', txt = 'gaussiano';
    case 'cir', txt = 'circular';
    otherwise, txt = tipo;
end
end


function [lambda, mu, varPred, rc, idx] = kriegeSistema(x, y, xq, yq, modelo, nMax)
%KRIEGESISTEMA Resuelve el sistema de Kriging Ordinario con los nMax
%vecinos mas cercanos al punto de consulta:
%   [Gamma 1; 1' 0] [lambda; mu] = [gamma0; 1],  sum(lambda) = 1
%   varianza de kriging = lambda' * gamma0 + mu
d = hypot(x - xq, y - yq);
[dOrd, orden] = sort(d);
k = min(nMax, numel(x));
idx = orden(1:k);
Xk = x(idx); Yk = y(idx);

Dij = hypot(Xk - Xk.', Yk - Yk.');
Gij = reshape(modeloVariograma(Dij(:), modelo.nugget, modelo.sill, modelo.rango, modelo.tipo), k, k);
A = [Gij, ones(k, 1); ones(1, k), 0];
g0 = modeloVariograma(dOrd(1:k), modelo.nugget, modelo.sill, modelo.rango, modelo.tipo);
b = [g0; 1];

rc = rcond(A);
estado = warning('off', 'MATLAB:nearlySingularMatrix');
estado2 = warning('off', 'MATLAB:singularMatrix');
sol = A \ b;
warning(estado); warning(estado2);

lambda = sol(1:k); mu = sol(end);
varPred = sum(lambda .* g0) + mu;
end


% ============================================================
% COMPARACION DE METODOS Y VALIDACION CRUZADA
% ============================================================

function C = compararMetodos(M, puntos, par, metodoRef)
%COMPARARMETODOS Evalua los seis metodos en cada punto y calcula la
%diferencia de cada uno respecto al metodo de referencia (m, cm y mm).
metodos = listaMetodos();
iRef = find(strcmp(metodos, metodoRef), 1);
nP = height(puntos); nM = numel(metodos);
Z = nan(nP, nM); W = nan(nP, nM); estado = strings(nP, nM);
for k = 1:nP
    for m = 1:nM
        r = interpolarPunto(M, puntos.Lat(k), puntos.Lon(k), metodos{m}, par, false);
        Z(k, m) = r.zeta; W(k, m) = r.sumaPesos;
        if r.aplicable
            estado(k, m) = "Calculado";
        else
            estado(k, m) = string(r.motivo);
        end
    end
end
D = Z - Z(:, iRef);

punto = repelem(puntos.Punto, nM, 1);
metodo = repmat(string(metodos(:)), nP, 1);
Zl = reshape(Z.', [], 1); Dl = reshape(D.', [], 1); Wl = reshape(W.', [], 1);
HN = repelem(puntos.h, nM, 1) - Zl;
C = struct();
C.tabla = table(punto, metodo, Zl, HN, Dl, 100*Dl, 1000*Dl, Wl, reshape(estado.', [], 1), ...
    'VariableNames', {'Punto', 'Metodo', 'Zeta_m', 'HN_m', 'Dif_m', 'Dif_cm', 'Dif_mm', ...
    'SumaPesos', 'Estado'});
C.Z = Z; C.D = D; C.metodos = metodos; C.referencia = metodoRef; C.nombres = puntos.Punto;
C.rangoMm = 1000 * (max(Z, [], 2) - min(Z, [], 2));
Dsin = abs(D); Dsin(isnan(Dsin)) = -Inf;
[mx, im] = max(Dsin, [], 2);
C.maxAbsMm = 1000 * mx;
C.metodoMaxAbs = string(metodos(im)).';
C.metodoMaxAbs = C.metodoMaxAbs(:);
end


function cv = validacionDecimacion(M, par, nControl, semilla, fnProgreso, metodos)
%VALIDACIONDECIMACION Validacion cruzada adaptada a mallas regulares.
%   El leave-one-out clasico no sirve para metodos exactos de malla (ver el
%   caso bilineal en completo6.m: al retirar el nodo, el punto consultado
%   cae justo sobre el hueco). Aqui se construye una malla decimada con los
%   nodos de indice impar (paso 2*dLat = 10') y se predicen, con cada metodo,
%   nodos retirados elegidos al azar con semilla fija (centros de celda y
%   puntos medios de aristas de la malla decimada). El error es
%   zeta_predicha - zeta_archivo. Como el paso de prueba es el doble del
%   real, los errores son una cota conservadora del error al paso de 5'.
if nargin < 5
    fnProgreso = [];
end
Mc = struct('lat', M.lat(1:2:end), 'lon', M.lon(1:2:end), 'Z', M.Z(1:2:end, 1:2:end), ...
    'dLat', 2*M.dLat, 'dLon', 2*M.dLon, 'nombre', [M.nombre ' (decimada)']);
[II, JJ] = ndgrid(1:numel(M.lat), 1:numel(M.lon));
cand = (mod(II - 1, 2) == 1 | mod(JJ - 1, 2) == 1) & ...
    M.lat(II) <= Mc.lat(end) & M.lon(JJ) <= Mc.lon(end) & ~isnan(M.Z);
idxCand = find(cand);
s = RandStream('mt19937ar', 'Seed', semilla);
nC = min(nControl, numel(idxCand));
sel = idxCand(randperm(s, numel(idxCand), nC));
latC = M.lat(II(sel)); lonC = M.lon(JJ(sel)); zC = M.Z(sel);

if nargin < 6 || isempty(metodos)
    metodos = listaMetodos();
end
nM = numel(metodos);
E = nan(nC, nM);
for m = 1:nM
    for k = 1:nC
        r = interpolarPunto(Mc, latC(k), lonC(k), metodos{m}, par, false);
        E(k, m) = r.zeta - zC(k);
        if ~isempty(fnProgreso) && mod(k, 25) == 0
            fnProgreso(((m - 1)*nC + k) / (nM*nC), sprintf('%s: %d/%d', metodos{m}, k, nC));
        end
    end
end

Emm = 1000 * E;
nEval = sum(~isnan(Emm), 1).';
cv = struct();
cv.tabla = table(string(metodos(:)), nEval, nC - nEval, mean(Emm, 1, 'omitnan').', ...
    sqrt(mean(Emm.^2, 1, 'omitnan')).', mean(abs(Emm), 1, 'omitnan').', ...
    percentilAbs(Emm, 95).', max(abs(Emm), [], 1).', ...
    'VariableNames', {'Metodo', 'Evaluados', 'NoAplicables', 'Sesgo_mm', 'RMSE_mm', ...
    'MAE_mm', 'P95_mm', 'MaxAbs_mm'});
cv.erroresMm = Emm; cv.lat = latC; cv.lon = lonC; cv.nControl = nC; cv.semilla = semilla;
cv.pasoPruebaMin = Mc.dLat * 60; cv.metodos = metodos;
end


function q = percentilAbs(E, pct)
%PERCENTILABS Percentil de |E| por columna (sin Statistics Toolbox;
%criterio del rango mas cercano).
q = nan(1, size(E, 2));
for c = 1:size(E, 2)
    v = sort(abs(E(~isnan(E(:, c)), c)));
    if ~isempty(v)
        q(c) = v(max(1, ceil(pct/100 * numel(v))));
    end
end
end


% ============================================================
% COMPARACION ENTRE MODELOS (QGeoidCOL2023, GeoCol2004, EGM2008)
% ============================================================

function [Tout, est] = compararModelos(Tin)
%COMPARARMODELOS Diferencias QGeoidCOL2023 - GeoCol2004 y QGeoidCOL2023 -
%EGM2008 (m, cm, mm) y alturas derivadas de cada modelo. Si h esta
%disponible: HN = h - zeta; H(GeoCol) = h - N; H(EGM) = h - N_EGM. El efecto
%del modelo sobre la altura es H(otro) - HN = zeta - N_otro = diferencia.
%Si Tin.H_nivelada (altura de la nivelacion geometrica) esta disponible, se
%calcula ademas el control altimetrico: H(modelo) - H_nivelada, en cm.
dGC = Tin.ZetaQ - Tin.N_GeoCol2004;
dEGM = Tin.ZetaQ - Tin.EGM2008;
HN = Tin.h - Tin.ZetaQ; HGC = Tin.h - Tin.N_GeoCol2004; HEGM = Tin.h - Tin.EGM2008;
Hn = Tin.H_nivelada;
Tout = table(Tin.Punto, Tin.Lat, Tin.Lon, Tin.h, Tin.ZetaQ, Tin.N_GeoCol2004, Tin.EGM2008, ...
    dGC, 100*dGC, 1000*dGC, dEGM, 100*dEGM, 1000*dEGM, HN, HGC, HEGM, ...
    Hn, 100*(HN - Hn), 100*(HGC - Hn), 100*(HEGM - Hn), ...
    'VariableNames', {'Punto', 'Lat', 'Lon', 'h', 'ZetaQ', 'N_GeoCol2004', 'EGM2008', ...
    'Q_GC_m', 'Q_GC_cm', 'Q_GC_mm', 'Q_EGM_m', 'Q_EGM_cm', 'Q_EGM_mm', ...
    'HN_m', 'H_GeoCol_m', 'H_EGM_m', 'H_niv_m', 'HN_menos_Hniv_cm', ...
    'HGC_menos_Hniv_cm', 'HEGM_menos_Hniv_cm'});
est = struct('GC', estadisticasDif(dGC), 'EGM', estadisticasDif(dEGM), ...
    'NivQ', estadisticasDif(HN - Hn), 'NivGC', estadisticasDif(HGC - Hn), ...
    'NivEGM', estadisticasDif(HEGM - Hn));
end


function s = estadisticasDif(d)
d = d(~isnan(d));
s = struct('n', numel(d), 'media', NaN, 'desv', NaN, 'rms', NaN, 'min', NaN, 'max', NaN);
if ~isempty(d)
    s.media = mean(d); s.rms = sqrt(mean(d.^2)); s.min = min(d); s.max = max(d);
    if numel(d) > 1
        s.desv = std(d);
    end
end
end


function lin = lineasCompatibilidad(egmEsZeta)
%LINEASCOMPATIBILIDAD Nota de compatibilidad conceptual de las magnitudes
%comparadas (requisito del enunciado antes de calcular diferencias).
lin = {'COMPARACION ENTRE MODELOS -- compatibilidad de magnitudes'; ...
    '  QGeoidCOL2023: anomalia de altura zeta (cuasigeoide, GRS80, IHRS, marea cero) -> HN = h - zeta.'; ...
    '  GeoCol2004: ondulacion geoidal N (geoide, GRS80, malla de 2'') -> H = h - N.'; ...
    '    Compatibilidad PARCIAL: zeta - N incluye la separacion geoide-cuasigeoide'; ...
    '    (N - zeta ~ Delta_g_Bouguer * H / gamma_medio; cm a dm en zonas de montana),'; ...
    '    ademas de las diferencias propias de los modelos (datos, metodo, datum vertical).'};
if egmEsZeta
    lin = [lin; {'  EGM2008 como anomalia de altura zeta: magnitud COMPATIBLE con QGeoidCOL2023.'; ...
        '    Verificar sistema de mareas y termino de grado cero de la herramienta usada.'}];
else
    lin = [lin; {'  EGM2008 como ondulacion geoidal N: compatibilidad PARCIAL (igual que GeoCol2004).'; ...
        '    Verificar sistema de mareas y termino de grado cero de la herramienta usada.'}];
end
lin = [lin; {'  QGeoidCOL2023 se adopta como referencia de comparacion, no como valor verdadero.'}];
end


function lin = lineasEstadisticasModelos(est)
lin = {'  Estadisticas de las diferencias (cm):'};
nombres = {'GC', 'QGeoid - GeoCol2004'; 'EGM', 'QGeoid - EGM2008'};
for k = 1:2
    s = est.(nombres{k, 1});
    if s.n == 0
        lin{end+1, 1} = sprintf('    %-20s sin valores externos', nombres{k, 2}); %#ok<AGROW>
    else
        lin{end+1, 1} = sprintf(['    %-20s n=%d  media=%+.1f  desv=%.1f  RMS=%.1f  ' ...
            'min=%+.1f  max=%+.1f'], nombres{k, 2}, s.n, 100*s.media, 100*s.desv, ...
            100*s.rms, 100*s.min, 100*s.max); %#ok<AGROW>
    end
end
nivel = {'NivQ', 'HN(QGeoid) - H nivelada'; 'NivGC', 'H(GeoCol2004) - H niv.'; ...
    'NivEGM', 'H(EGM2008) - H niv.'};
if est.NivQ.n > 0
    lin{end+1, 1} = '  Control altimetrico frente a la nivelacion (cm):';
    for k = 1:3
        s = est.(nivel{k, 1});
        lin{end+1, 1} = sprintf(['    %-24s n=%d  media=%+.1f  desv=%.1f  RMS=%.1f  ' ...
            'min=%+.1f  max=%+.1f'], nivel{k, 2}, s.n, 100*s.media, 100*s.desv, ...
            100*s.rms, 100*s.min, 100*s.max); %#ok<AGROW>
    end
end
end


function [Tin, nAsignados] = importarValoresExternos(ruta, Tin)
%IMPORTARVALORESEXTERNOS Lee un CSV/TXT con columnas cuyos nombres
%contengan "geocol" y "egm" (y opcionalmente "punto" o "nombre" y una
%columna con "nivelada", la altura de la nivelacion geometrica). Si hay columna de punto se asigna
%por nombre; si no, por orden de fila.
Text = readtable(ruta, 'TextType', 'string', 'VariableNamingRule', 'preserve', ...
    'VariableNamesLine', 1);
vars = lower(string(Text.Properties.VariableNames));
iGC = find(contains(vars, 'geocol'), 1);
iE = find(contains(vars, 'egm'), 1);
iP = find(vars == "punto" | contains(vars, 'nombre'), 1);
iH = find(contains(vars, 'nivelada'), 1);
if isempty(iGC) || isempty(iE)
    error('importarValoresExternos:formato', ['El archivo debe tener columnas cuyo nombre ' ...
        'contenga "GeoCol" y "EGM" (p.ej. Punto,GeoCol2004,EGM2008).']);
end
nAsignados = 0;
for k = 1:height(Tin)
    if ~isempty(iP)
        m = find(strtrim(string(Text{:, iP})) == Tin.Punto(k), 1);
    elseif k <= height(Text)
        m = k;
    else
        m = [];
    end
    if ~isempty(m)
        Tin.N_GeoCol2004(k) = double(Text{m, iGC});
        Tin.EGM2008(k) = double(Text{m, iE});
        if ~isempty(iH)
            Tin.H_nivelada(k) = double(Text{m, iH});
        end
        nAsignados = nAsignados + 1;
    end
end
end


% ============================================================
% UTILIDADES DE LA INTERFAZ (LOG, PARSEO DE PUNTOS, TABLAS, TEXTO)
% ============================================================

function app = agregarLog(app, texto)
%AGREGARLOG Anexa una linea (o varias) al panel de bitacora de la app.
if ischar(texto) || isstring(texto)
    texto = cellstr(texto);
end
% [MODIFICADO] formato de la bitacora: se descartan las lineas de guiones,
% se antepone la hora al bloque, se sangran las lineas siguientes y se deja
% una linea en blanco al final; luego se baja al final del cuadro.
texto = texto(:);
esSep = cellfun(@(s) ischar(s) && ~isempty(s) && all(s == '-'), texto);
texto = texto(~esSep);
if ~isempty(texto)
    sello = ['[' char(datetime('now', 'Format', 'HH:mm:ss')) ']  '];
    sangria = repmat(' ', 1, numel(sello));
    texto{1} = [sello texto{1}];
    for k = 2:numel(texto)
        texto{k} = [sangria texto{k}];
    end
    texto{end+1, 1} = '';
end
app.txtLog.Value = [app.txtLog.Value; texto];
try
    scroll(app.txtLog, 'bottom');
catch
end
drawnow limitrate;
end


function actualizarProgreso(dlg, fraccion, msg)
dlg.Value = min(max(fraccion, 0), 1);
dlg.Message = msg;
end


function T = tablaControlCalidad(rep, verif)
%TABLACONTROLCALIDAD Arma la tabla (Grupo, Control, Resultado, Estado,
%Criterio) del control de calidad de la malla y la autoverificacion.
E = {'REVISAR', 'OK'};
SN = {'NO', 'SI'};
tol = 1e-9;
L = 'Lectura del archivo'; S = 'Estructura de la malla'; C = 'Cobertura'; Z = 'Anomalia zeta';
T = { ...
    L, 'Lineas con contenido', sprintf('%d', rep.nLineas), '-', 'Cada linea debe traer latitud, longitud y zeta'; ...
    L, 'Lineas sin 3 numeros', sprintf('%d', rep.nNoNumericas), E{(rep.nNoNumericas == 0) + 1}, 'Debe ser 0'; ...
    L, 'Valores NaN/Inf', sprintf('%d', rep.nNoFinitos), E{(rep.nNoFinitos == 0) + 1}, 'Debe ser 0'; ...
    L, 'Coordenadas fuera de rango', sprintf('%d', rep.nFueraRango), E{(rep.nFueraRango == 0) + 1}, 'Latitud entre -90 y 90, longitud entre -180 y 180. Debe ser 0'; ...
    L, 'Nodos duplicados', sprintf('%d', rep.nDuplicados), E{(rep.nDuplicados == 0) + 1}, 'Misma posicion repetida. Debe ser 0'; ...
    L, 'Duplicados con zeta distinta', sprintf('%d', rep.nDuplicadosConflicto), E{(rep.nDuplicadosConflicto == 0) + 1}, 'Debe ser 0'; ...
    S, 'Dimensiones (lat x lon)', sprintf('%d x %d', rep.nLat, rep.nLon), '-', 'Numero de latitudes por numero de longitudes'; ...
    S, 'Nodos esperados', sprintf('%d', rep.nEsperados), '-', 'Latitudes x longitudes'; ...
    S, 'Nodos presentes / faltantes', sprintf('%d / %d', rep.nPresentes, rep.nFaltantes), E{(rep.nFaltantes == 0) + 1}, 'Malla rectangular completa: no debe faltar ningun nodo'; ...
    S, 'Malla regular (paso fijo)', SN{rep.esRegular + 1}, E{rep.esRegular + 1}, 'El paso debe ser constante en latitud y en longitud'; ...
    S, 'Paso en latitud', sprintf('%.4f min (%.10f grados)', rep.dLat*60, rep.dLat), '-', sprintf('Variacion maxima del paso: %.1e grados', rep.desvPasoLat); ...
    S, 'Paso en longitud', sprintf('%.4f min (%.10f grados)', rep.dLon*60, rep.dLon), '-', sprintf('Variacion maxima del paso: %.1e grados', rep.desvPasoLon); ...
    S, 'Desviacion vs malla ideal', sprintf('%.1e grados', rep.desvNominal), '-', 'Distancia maxima entre cada nodo del archivo y su posicion teorica'; ...
    S, 'Resolucion', sprintf('%.3f km (N-S) x %.3f km (E-O)', rep.resLatKm, rep.resLonKm), '-', 'Separacion entre nodos en el centro de la malla'; ...
    S, 'Orden del archivo', rep.orden, '-', 'Solo informativo: la malla se reconstruye por indices'; ...
    C, 'Latitud', sprintf('de %.6f a %.6f grados', rep.latLim(1), rep.latLim(2)), '-', 'Fuera de este rango no se interpola'; ...
    C, 'Longitud', sprintf('de %.6f a %.6f grados', rep.lonLim(1), rep.lonLim(2)), '-', 'Fuera de este rango no se interpola'; ...
    C, 'Extension aproximada', sprintf('%.0f km (N-S) x %.0f km (E-O)', rep.extensionNSkm, rep.extensionEOkm), '-', 'Tamano de la zona cubierta'; ...
    Z, 'zeta minimo', sprintf('%.3f m (lat %.4f, lon %.4f)', rep.zMin, rep.posMin), '-', 'Valor mas bajo de la malla y su posicion'; ...
    Z, 'zeta maximo', sprintf('%.3f m (lat %.4f, lon %.4f)', rep.zMax, rep.posMax), '-', 'Valor mas alto de la malla y su posicion'; ...
    Z, 'Media / desv. estandar', sprintf('%.3f m / %.3f m', rep.zMedia, rep.zDesv), '-', 'Estadistica de todos los nodos'; ...
    Z, '|zeta| mayor a 150 m', sprintf('%d', rep.nZetaImplausible), E{(rep.nZetaImplausible == 0) + 1}, 'Valores no plausibles. Debe ser 0'};
for m = 1:numel(verif.metodo)
    T(end+1, :) = {'Autoverificacion', ['zeta ' verif.metodo{m} ' vs interp2'], ...
        sprintf('%.2e m', verif.maxDif(m)), E{(verif.maxDif(m) < tol) + 1}, ...
        sprintf('Diferencia maxima en %d puntos; debe ser menor a 1e-9 m', verif.n)}; %#ok<AGROW>
end
for m = 1:numel(verif.metodosPesos)
    T(end+1, :) = {'Autoverificacion', ['Suma de pesos: ' verif.metodosPesos{m}], ...
        sprintf('%.1e', verif.maxErrSumaPesos(m)), E{(verif.maxErrSumaPesos(m) < tol) + 1}, ...
        sprintf('max|1 - suma(w)| en %d puntos; debe ser menor a 1e-9', verif.nPesos)}; %#ok<AGROW>
end
end


function lin = infoBase()
lin = { ...
    'MODELO: QGeoidCOL2023, anomalia de altura zeta (m) en una malla regular.'; ...
    'ARCHIVO: .txt con 3 columnas: latitud, longitud (grados decimales) y zeta (m).'; ...
    'USO: altura normal HN = h - zeta, con h = altura elipsoidal GNSS.'; ...
    ' '; ...
    'METODOS DE INTERPOLACION:'; ...
    '1. Vecino mas cercano: usa el zeta del nodo mas proximo.'; ...
    '2. Lineal por triangulos: plano sobre el triangulo (3 nodos) que contiene el punto.'; ...
    '3. Bilineal: promedio ponderado de los 4 nodos de la celda.'; ...
    '4. Bicubica: convolucion cubica (Keys) con 16 nodos (4x4).'; ...
    '5. IDW: promedio ponderado de los k nodos mas cercanos, peso = 1/d^p.'; ...
    '6. Kriging ordinario: pesos obtenidos del variograma ajustado, con n vecinos.'; ...
    ' '; ...
    'Los metodos 1-4 operan sobre los indices de la malla (lat/lon); IDW y kriging usan distancias en Origen-Nacional MAGNA-SIRGAS (EPSG:9377, km).'; ...
    'No se extrapola fuera de la malla. Presentacion: zeta y HN redondeados a 1 mm.'};
end


function puntos = parsearPuntos(lineasCelda)
%PARSEARPUNTOS Convierte el contenido del cuadro de texto en una tabla
%(Punto, Lat, Lon, h). Cada entrada: "[Nombre:] lat, lon [, h]", una por
%linea o separadas por ";". Lineas que empiezan por "%" o "#" se ignoran.
%Los numeros se extraen por expresion regular, asi que se toleran comas,
%espacios o ambos como separador (con punto decimal).
lineas = cellstr(lineasCelda);
entradas = {};
for i = 1:numel(lineas)
    s = strtrim(lineas{i});
    if isempty(s) || s(1) == '%' || s(1) == '#'
        continue
    end
    entradas = [entradas, strsplit(s, ';')]; %#ok<AGROW>
end

nombres = strings(0, 1); lat = zeros(0, 1); lon = zeros(0, 1); h = zeros(0, 1);
for i = 1:numel(entradas)
    s = strtrim(entradas{i});
    nombre = "";
    c = strfind(s, ':');
    if ~isempty(c)
        nombre = string(strtrim(s(1:c(1)-1)));
        s = s(c(1)+1:end);
    end
    nums = regexp(s, '[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?', 'match');
    if numel(nums) < 2
        continue
    end
    v = str2double(nums);
    if any(isnan(v(1:2))) || abs(v(1)) > 90 || abs(v(2)) > 180
        continue
    end
    if nombre == ""
        nombre = "P" + (numel(lat) + 1);
    end
    nombres(end+1, 1) = nombre; %#ok<AGROW>
    lat(end+1, 1) = v(1); %#ok<AGROW>
    lon(end+1, 1) = v(2); %#ok<AGROW>
    if numel(v) >= 3
        h(end+1, 1) = v(3); %#ok<AGROW>
    else
        h(end+1, 1) = NaN; %#ok<AGROW>
    end
end
puntos = table(nombres, lat, lon, h, 'VariableNames', {'Punto', 'Lat', 'Lon', 'h'});
end


function lineas = puntosControlPredeterminados()
%PUNTOSCONTROLPREDETERMINADOS Cinco puntos de control P1-P5 del archivo
%"5 PUNTOS CALCULADOS COLOMBIA.xlsx" (nomenclatura, latitud y longitud en
%grados decimales, altura elipsoidal h en m): Guajira, Cordoba, Norte de
%Santander, Bogota y Arauca. Todos caen dentro de la cobertura del modelo.
% [MODIFICADO] mismos datos, presentados en columnas alineadas (el formato
% "Nombre: lat, lon, h" y las lineas de comentario "%" siguen siendo validos).
nom = {'P1 32-TN-7'; 'P2 23555011'; 'P3 GPS-S-T-144'; 'P4 (EA)-MACARENA-CD-551-1'; 'P5 81794007'};
dat = [11.4519608895, -72.5744326992,   35.134; ...
        8.431102848,  -75.553010881,    98.266; ...
        7.311055511,  -72.8720174438, 3228.966; ...
        4.6149699549, -74.0628738019, 2725.381; ...
        6.2760436758, -71.6514462724,  185.792];
lineas = cell(numel(nom) + 3, 1);
lineas{1} = '% PUNTOS DE CONTROL (lat, lon en grados; h = alt. elipsoidal, m)';
lineas{2} = sprintf('%% %-24s %13s  %14s  %9s', 'Nombre', 'Lat', 'Lon', 'h (m)');
lineas{3} = ['% ' repmat('-', 1, 65)];
for k = 1:numel(nom)
    lineas{k + 3} = sprintf('%-26s %13.10f, %14.10f, %9.3f', [nom{k} ':'], dat(k, 1), dat(k, 2), dat(k, 3));
end
end


function T = tablaResultadosPuntos(puntos, resultados)
%TABLARESULTADOSPUNTOS Tabla resumen (precision completa) de los puntos
%consultados con el metodo activo.
n = height(puntos);
zeta = nan(n, 1); nNod = zeros(n, 1); sumW = nan(n, 1); sigma = nan(n, 1);
estado = strings(n, 1); obs = strings(n, 1);
for k = 1:n
    r = resultados{k};
    zeta(k) = r.zeta; sumW(k) = r.sumaPesos; sigma(k) = r.sigma;
    if r.aplicable
        estado(k) = "Calculado"; nNod(k) = height(r.nodos); obs(k) = string(r.info);
    else
        estado(k) = "Sin valor"; obs(k) = string(r.motivo);
    end
end
T = table(puntos.Punto, puntos.Lat, puntos.Lon, estado, zeta, puntos.h, puntos.h - zeta, ...
    nNod, sumW, sigma, obs, 'VariableNames', {'Punto', 'Lat', 'Lon', 'Estado', 'Zeta_m', ...
    'h_m', 'HN_m', 'NNodos', 'SumaPesos', 'SigmaKO_m', 'Observacion'});
end


function T = redondearParaPresentacion(T)
%REDONDEARPARAPRESENTACION Aplica el criterio de redondeo declarado (solo
%para mostrar): coordenadas a 1e-6 grados, alturas/zeta/diferencias en m a
%1 mm, en cm a 0.1 mm, en mm a 0.01 mm, sumas de pesos a 1e-12.
for v = string(T.Properties.VariableNames)
    if ~isnumeric(T.(v))
        continue
    end
    if any(v == ["Lat", "Lon"])
        T.(v) = round(T.(v), 6);
    elseif endsWith(v, "_cm")
        T.(v) = round(T.(v), 2);
    elseif endsWith(v, "_mm")
        T.(v) = round(T.(v), 2);
    elseif v == "SumaPesos"
        T.(v) = round(T.(v), 12);
    elseif v == "SigmaKO_m"
        T.(v) = round(T.(v), 4);
    elseif v == "NNodos"
        % entero
    else
        T.(v) = round(T.(v), 3);
    end
end
end


function lin = lineasResumenPunto(p, r, par)
lin = { ...
    sprintf('Punto:   %s', p.Punto); ...
    sprintf('Lat/Lon: %.6f, %.6f (grados)', p.Lat, p.Lon); ...
    sprintf('Metodo:  %s (%s)', r.metodo, textoParametros(par, r.metodo)); ...
    ''};
if ~r.aplicable
    lin = [lin; {'SIN VALOR:'}; {r.motivo}];
    return
end
errSuma = abs(1 - r.sumaPesos);
if errSuma < 1e-9
    okSuma = 'OK';
else
    okSuma = 'REVISAR';
end
lin = [lin; { ...
    sprintf('zeta interpolada  = %.3f m', r.zeta); ...
    sprintf('  (sin redondear  = %.12f m)', r.zeta); ...
    sprintf('h suministrada    = %s', textoNum(p.h, '%.3f m')); ...
    sprintf('HN = h - zeta     = %s', textoNum(p.h - r.zeta, '%.3f m')); ...
    ''; ...
    sprintf('Nodos usados: %d', height(r.nodos)); ...
    sprintf('Suma de pesos     = %.15f', r.sumaPesos); ...
    sprintf('|1 - suma|        = %.1e  -> %s', errSuma, okSuma); ...
    sprintf('Suma contribuc.   = %.12f m', sum(r.nodos.Contrib)); ...
    sprintf('Pesos negativos   = %d', sum(r.nodos.Peso < 0))}];
if ~isnan(r.sigma)
    lin{end+1} = sprintf('sigma kriging     = %.4f m', r.sigma);
end
lin = [lin; {''; ['Info: ' r.info]; ''; ...
    'Redondeo de presentacion: zeta y HN a 0.001 m;'; ...
    'pesos a 1e-6; contribuciones a 1e-4 m. Calculo en'; ...
    'doble precision (CSV exportados sin redondeo).'}];
end


function txt = textoNum(x, fmt)
if isnan(x)
    txt = 'no suministrada';
else
    txt = sprintf(fmt, x);
end
end


function txt = textoLineas(v)
if isempty(v)
    txt = '';
else
    txt = ['(primeras lineas: ' strjoin(arrayfun(@num2str, v, 'UniformOutput', false), ', ') ')'];
end
end


function txt = nombreArchivo(s)
txt = lower(regexprep(char(s), '[^A-Za-z0-9]+', '_'));
end


% ============================================================
% GRAFICOS
% ============================================================

function graficarMallaGeneral(ax, M, puntos)
%GRAFICARMALLAGENERAL Vista general de la malla QGeoidCOL2023 (zeta en
%color, isolineas cada 5 m, limite de cobertura) con los puntos consultados.
cla(ax, 'reset');
imagesc(ax, M.lon, M.lat, M.Z, 'AlphaData', ~isnan(M.Z));
ax.YDir = 'normal';
colormap(ax, turbo(256));
cb = colorbar(ax);
cb.Label.String = '\zeta (m)';
hold(ax, 'on');
niveles = 5*floor(min(M.Z(:))/5):5:5*ceil(max(M.Z(:))/5);
contour(ax, M.lon, M.lat, M.Z, niveles, 'k', 'LineWidth', 0.4, 'ShowText', 'on', ...
    'LabelSpacing', 500);
plot(ax, M.lon([1 end end 1 1]), M.lat([1 1 end end 1]), 'k-', 'LineWidth', 1.2);
if height(puntos) > 0
    plot(ax, puntos.Lon, puntos.Lat, 'p', 'MarkerSize', 13, 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    text(ax, puntos.Lon, puntos.Lat, "  " + puntos.Punto, 'FontWeight', 'bold', ...
        'FontSize', 9, 'Color', 'k', 'BackgroundColor', [1 1 1 0.6], 'Margin', 0.5);
end
hold(ax, 'off');
daspect(ax, [1/cosd(mean(M.lat([1 end]))) 1 1]);
xlim(ax, M.lon([1 end]).' + [-0.2 0.2]); ylim(ax, M.lat([1 end]).' + [-0.2 0.2]);
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'Longitud (grados)'); ylabel(ax, 'Latitud (grados)');
title(ax, sprintf('%s: %d x %d nodos, paso %.0f''', M.nombre, numel(M.lat), numel(M.lon), M.dLat*60), ...
    'Interpreter', 'none');
end


function graficarHistogramaZeta(ax, M)
cla(ax, 'reset');
histogram(ax, M.Z(:), 60, 'FaceColor', [0.15 0.40 0.75], 'EdgeColor', 'none');
grid(ax, 'on');
xlabel(ax, '\zeta (m)'); ylabel(ax, 'Numero de nodos');
title(ax, sprintf('Histograma de \\zeta (media %.2f m, desv. %.2f m)', ...
    mean(M.Z(:), 'omitnan'), std(M.Z(:), 'omitnan')));
end


function campo = campoLocal(M, r, par)
%CAMPOLOCAL Evalua el metodo activo en una rejilla fina alrededor del
%punto (+-2.5 nodos) para la vista local. En kriging se reutiliza el
%variograma ajustado para el punto consultado (un solo ajuste).
campo = [];
if ~dentroDeCobertura(M, r.lat, r.lon)
    return
end
semi = 2.5;
latLim = [max(M.lat(1), r.lat - semi*M.dLat), min(M.lat(end), r.lat + semi*M.dLat)];
lonLim = [max(M.lon(1), r.lon - semi*M.dLon), min(M.lon(end), r.lon + semi*M.dLon)];
nG = 35;
latV = linspace(latLim(1), latLim(2), nG);
lonV = linspace(lonLim(1), lonLim(2), nG);
if strcmp(r.metodo, 'Kriging ordinario')
    if isempty(r.modeloKrig)
        return
    end
    par.modeloKrigFijo = r.modeloKrig;
end
Zc = nan(nG, nG);
for a = 1:nG
    for b = 1:nG
        Zc(a, b) = interpolarPunto(M, latV(a), lonV(b), r.metodo, par, false).zeta;
    end
end
campo = struct('lat', latV, 'lon', lonV, 'Z', Zc, 'latLim', latLim, 'lonLim', lonLim);
end


function graficarVistaLocal(ax, M, r, campo, nombrePunto)
%GRAFICARVISTALOCAL Campo interpolado alrededor del punto, nodos de la
%malla, nodos usados (tamano ~ |peso|, rojo w>0, azul w<0, etiqueta con w)
%y geometria de soporte del metodo.
cla(ax, 'reset');
hold(ax, 'on');
leyH = gobjects(0); leyT = {};
if ~isempty(campo) && any(~isnan(campo.Z(:)))
    contourf(ax, campo.lon, campo.lat, campo.Z, 24, 'LineColor', 'none');
    colormap(ax, turbo(256));
    cb = colorbar(ax); cb.Label.String = '\zeta interpolada (m)';
    latLim = campo.latLim; lonLim = campo.lonLim;
else
    latLim = r.lat + [-2.5 2.5]*M.dLat; lonLim = r.lon + [-2.5 2.5]*M.dLon;
end

% Nodos de la malla en la ventana visible
ii = find(M.lat >= latLim(1) - 1.5*M.dLat & M.lat <= latLim(2) + 1.5*M.dLat);
jj = find(M.lon >= lonLim(1) - 1.5*M.dLon & M.lon <= lonLim(2) + 1.5*M.dLon);
[LONn, LATn] = meshgrid(M.lon(jj), M.lat(ii));
leyH(end+1) = plot(ax, LONn(:), LATn(:), 'k.', 'MarkerSize', 8); leyT{end+1} = 'Nodos de la malla';

% Geometria de soporte
g = r.geom;
if ~isempty(g.celdaLat)
    plot(ax, g.celdaLon, g.celdaLat, 'w--', 'LineWidth', 1.2);
end
switch g.tipo
    case 'vecino'
        leyH(end+1) = plot(ax, g.lon, g.lat, 'w-', 'LineWidth', 2); leyT{end+1} = 'Soporte';
    case 'poligono'
        if ~isempty(g.lat)
            leyH(end+1) = plot(ax, g.lon, g.lat, 'w-', 'LineWidth', 2.2); leyT{end+1} = 'Soporte';
        end
    case 'radial'
        if ~isempty(g.lat)
            xs = [repmat(r.lon, 1, numel(g.lon)); g.lon(:).'; nan(1, numel(g.lon))];
            ys = [repmat(r.lat, 1, numel(g.lat)); g.lat(:).'; nan(1, numel(g.lat))];
            leyH(end+1) = plot(ax, xs(:), ys(:), 'w-', 'LineWidth', 0.9); leyT{end+1} = 'Soporte';
        end
end

% Nodos usados con su peso
if r.aplicable
    N = r.nodos;
    tam = 30 + 320 * abs(N.Peso) / max(abs(N.Peso));
    col = repmat([0.85 0.10 0.10], height(N), 1);
    col(N.Peso < 0, :) = repmat([0.10 0.30 0.90], sum(N.Peso < 0), 1);
    leyH(end+1) = scatter(ax, N.Lon, N.Lat, tam, col, 'filled', 'MarkerEdgeColor', 'k');
    leyT{end+1} = 'Nodos usados (tamano ~ |w|; azul: w < 0)';
    if height(N) <= 25
        text(ax, N.Lon, N.Lat, compose('  %.3f', N.Peso), 'FontSize', 8, ...
            'Color', 'k', 'BackgroundColor', [1 1 1 0.65], 'Margin', 0.5);
    end
end
leyH(end+1) = plot(ax, r.lon, r.lat, 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.3);
leyT{end+1} = 'Punto consultado';
hold(ax, 'off');

xlim(ax, lonLim); ylim(ax, latLim);
daspect(ax, [1/cosd(r.lat) 1 1]);
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'Longitud (grados)'); ylabel(ax, 'Latitud (grados)');
if r.aplicable
    title(ax, sprintf('%s -- %s: \\zeta = %.3f m', nombrePunto, r.metodo, r.zeta));
else
    title(ax, sprintf('%s -- %s: sin valor', nombrePunto, r.metodo));
end
legend(ax, leyH, leyT, 'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 8);
end


function graficarComparacionMetodos(ax, C)
cla(ax, 'reset');
n = numel(C.nombres);
b = bar(ax, 1:n, 1000 * C.D, 'grouped');
cols = lines(numel(C.metodos));
for k = 1:numel(b)
    b(k).FaceColor = cols(k, :);
end
xticks(ax, 1:n); xticklabels(ax, cellstr(C.nombres));
xlabel(ax, 'Punto de consulta');
ylabel(ax, sprintf('Diferencia de \\zeta respecto a %s (mm)', C.referencia));
grid(ax, 'on');
lg = legend(ax, C.metodos, 'Location', 'eastoutside');
title(lg, 'Metodo de interpolacion');
title(ax, sprintf('Diferencia de cada metodo respecto a la referencia (%s)', C.referencia));
end


function graficarValidacionCruzada(ax, axBar, cv)
%GRAFICARVALIDACIONCRUZADA Diagrama de cajas de los errores (eje recortado
%al mayor P95 para que los metodos de menor error sigan siendo legibles) y
%barras de RMSE, MAE, P95 y maximo en escala logaritmica.
cla(ax, 'reset');
nC = size(cv.erroresMm, 1);
e = reshape(cv.erroresMm.', [], 1);
grupos = categorical(repmat(string(cv.metodos(:)), nC, 1), string(cv.metodos(:)));
ok = ~isnan(e);
boxchart(ax, grupos(ok), e(ok), 'MarkerStyle', '.');
yline(ax, 0, 'k-');
lim = max(cv.tabla.P95_mm);
ylim(ax, [-lim lim]);
grid(ax, 'on');
ylabel(ax, 'Error de prediccion (mm)');
title(ax, sprintf('Errores (%d nodos, malla a %.0f''), eje a \\pm%.0f mm', ...
    cv.nControl, cv.pasoPruebaMin, lim));

cla(axBar, 'reset');
nM = numel(cv.metodos);
b = bar(axBar, 1:nM, [cv.tabla.RMSE_mm, cv.tabla.MAE_mm, cv.tabla.P95_mm, cv.tabla.MaxAbs_mm], 'grouped');
cols = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.93 0.69 0.13; 0.64 0.08 0.18];
for k = 1:numel(b)
    b(k).FaceColor = cols(k, :);
end
axBar.YScale = 'log';
xticks(axBar, 1:nM); xticklabels(axBar, cv.metodos);
ylabel(axBar, 'mm (escala logaritmica)');
grid(axBar, 'on');
legend(axBar, {'RMSE', 'MAE', 'P95', 'Maximo'}, 'Location', 'northeast');
title(axBar, 'Indicadores de error por metodo');
end


function graficarDiferenciasModelos(axGC, axEGM, axBar, M, T)
%GRAFICARDIFERENCIASMODELOS Comportamiento espacial (puntos coloreados por
%la diferencia en cm sobre isolineas de zeta) y diagrama de barras.
datos = {axGC, T.Q_GC_m, 'QGeoidCOL2023 - GeoCol2004'; ...
         axEGM, T.Q_EGM_m, 'QGeoidCOL2023 - EGM2008'};
niveles = 5*floor(min(M.Z(:))/5):5:5*ceil(max(M.Z(:))/5);
for k = 1:2
    ax = datos{k, 1}; d = 100 * datos{k, 2};
    cla(ax, 'reset');
    contour(ax, M.lon, M.lat, M.Z, niveles, 'Color', [0.7 0.7 0.7]);
    hold(ax, 'on');
    ok = ~isnan(d);
    if any(ok)
        scatter(ax, T.Lon(ok), T.Lat(ok), 110, d(ok), 'filled', 'MarkerEdgeColor', 'k');
        colormap(ax, mapaDivergente(256));
        lim = max(abs(d(ok))); if lim == 0, lim = 1; end
        clim(ax, [-lim lim]);
        cb = colorbar(ax); cb.Label.String = 'Diferencia (cm)';
        % Solo la primera palabra del nombre (p.ej. "P1") para que la
        % etiqueta no se salga del mapa.
        idCorto = extractBefore(T.Punto(ok) + " ", " ");
        text(ax, T.Lon(ok), T.Lat(ok), compose('  %s: %+.1f cm', idCorto, d(ok)), ...
            'FontSize', 8, 'BackgroundColor', [1 1 1 0.7], 'Margin', 0.5);
    end
    plot(ax, T.Lon(~ok), T.Lat(~ok), 'kx', 'MarkerSize', 9);
    hold(ax, 'off');
    xlim(ax, M.lon([1 end])); ylim(ax, M.lat([1 end]));
    daspect(ax, [1/cosd(mean(M.lat([1 end]))) 1 1]);
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'Longitud (grados)'); ylabel(ax, 'Latitud (grados)');
    title(ax, datos{k, 3});
end

cla(axBar, 'reset');
n = height(T);
b = bar(axBar, 1:n, 100 * [T.Q_GC_m, T.Q_EGM_m], 'grouped');
b(1).FaceColor = [0.85 0.33 0.10]; b(2).FaceColor = [0.00 0.45 0.74];
xticks(axBar, 1:n); xticklabels(axBar, cellstr(T.Punto));
ylabel(axBar, 'Diferencia (cm)'); grid(axBar, 'on');
% Margen superior para que la leyenda no tape las barras.
vals = 100 * [T.Q_GC_m; T.Q_EGM_m]; vals = vals(~isnan(vals));
if ~isempty(vals)
    lo = min([vals; 0]); hi = max([vals; 0]); rng = max(hi - lo, 1);
    ylim(axBar, [lo - 0.05*rng, hi + 0.35*rng]);
end
legend(axBar, {'QGeoid - GeoCol2004', 'QGeoid - EGM2008'}, 'Location', 'northwest');
title(axBar, 'Diferencias entre modelos por punto (cm)');
end


function cmap = mapaDivergente(n)
%MAPADIVERGENTE Mapa de color azul-blanco-rojo para diferencias con signo.
m = floor(n/2);
azul = [0.13 0.40 0.67]; rojo = [0.70 0.09 0.17]; blanco = [1 1 1];
t1 = linspace(0, 1, m).'; t2 = linspace(0, 1, n - m).';
cmap = [azul + t1*(blanco - azul); blanco + t2*(rojo - blanco)];
end


% [NUEVO] ============================================================
% VARIOGRAMA Y COVARIOGRAMA USADOS POR EL KRIGING ORDINARIO
% ====================================================================

function graficarVariograma(axG, axC, tbl, r, par, nombrePunto)
%GRAFICARVARIOGRAMA Muestra, para el punto en detalle, el semivariograma
%experimental (puntos, tamano ~ numero de pares), el modelo ajustado por
%MCP, el covariograma derivado C(h) = (C0 + C1) - gamma(h) para h > 0, con
%C(0) = C0 + C1, y la tabla de parametros. Solo aplica con "Kriging
%ordinario"; con otro metodo se informa que no hay variograma.
cla(axG, 'reset'); cla(axC, 'reset');
if ~strcmp(r.metodo, 'Kriging ordinario') || isempty(r.modeloKrig)
    title(axG, 'Semivariograma ajustado: sin datos');
    title(axC, 'Covariograma derivado: sin datos');
    if strcmp(r.metodo, 'Kriging ordinario')
        msg = 'No se pudo ajustar el variograma para este punto (ver "Observacion" en la tabla de puntos).';
    else
        msg = ['El metodo activo (' r.metodo ') no usa variograma. ' ...
            'Seleccione "Kriging ordinario" y presione "Calcular".'];
    end
    tbl.Data = {'Estado', msg, ''};
    return
end

m = r.modeloKrig;
cTot = m.nugget + m.sill;                 % C(0) = C0 + C1
switch m.tipo                             % rango practico (95 % de la meseta)
    case 'exp', rPrac = 3*m.rango;
    case 'gau', rPrac = sqrt(3)*m.rango;
    otherwise,  rPrac = m.rango;          % esferico y circular: el rango es el practico
end
hh = linspace(0, m.lagMax, 200).';
gm = modeloVariograma(hh, m.nugget, m.sill, m.rango, m.tipo);
cm = cTot - gm;
tam = 25 + 175 * m.nPares(:) / max(m.nPares);
azul = [0.00 0.45 0.74];

% --- Semivariograma ---
hold(axG, 'on');
leyH = gobjects(0); leyT = {};
leyH(end+1) = scatter(axG, m.hExp, m.gammaExp, tam, azul, 'filled', 'MarkerEdgeColor', 'k');
leyT{end+1} = 'Experimental (tamano ~ no. de pares)';
leyH(end+1) = plot(axG, hh(2:end), gm(2:end), 'r-', 'LineWidth', 2);
leyT{end+1} = ['Modelo ' nombreModelo(m.tipo) ' (MCP)'];
leyH(end+1) = yline(axG, cTot, 'k--');
leyT{end+1} = 'Meseta total C0 + C1';
leyH(end+1) = yline(axG, m.nugget, ':', 'Color', [0.4 0.4 0.4]);
leyT{end+1} = 'Pepita C0';
hold(axG, 'off');
xlim(axG, [0 m.lagMax]);
grid(axG, 'on'); box(axG, 'on');
xlabel(axG, 'Distancia h (km)'); ylabel(axG, 'Semivarianza \gamma(h) (m^2)');
title(axG, sprintf('%s -- semivariograma %s', nombrePunto, nombreModelo(m.tipo)), 'Interpreter', 'none');
legend(axG, leyH, leyT, 'Location', 'southoutside', 'FontSize', 8);

% --- Covariograma derivado ---
hold(axC, 'on');
leyH = gobjects(0); leyT = {};
leyH(end+1) = scatter(axC, m.hExp, cTot - m.gammaExp, tam, azul, 'filled', 'MarkerEdgeColor', 'k');
leyT{end+1} = 'Experimental: (C0 + C1) - \gamma_{exp}';
leyH(end+1) = plot(axC, hh(2:end), cm(2:end), 'r-', 'LineWidth', 2);
leyT{end+1} = 'Modelo: C(h) = (C0 + C1) - \gamma(h)';
leyH(end+1) = plot(axC, 0, cTot, 's', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
leyT{end+1} = 'C(0) = C0 + C1';
hold(axC, 'off');
xlim(axC, [0 m.lagMax]);
grid(axC, 'on'); box(axC, 'on');
xlabel(axC, 'Distancia h (km)'); ylabel(axC, 'Covarianza C(h) (m^2)');
title(axC, sprintf('%s -- covariograma derivado', nombrePunto), 'Interpreter', 'none');
legend(axC, leyH, leyT, 'Location', 'southoutside', 'FontSize', 8);

% --- Tabla de parametros ---
semiEfectivo = max(par.ventanaKrig, ceil(sqrt(par.nKrig)/2) + 1);
if r.aplicable
    nVec = height(r.nodos);
else
    nVec = par.nKrig;
end
tbl.Data = { ...
    'Punto en detalle', char(nombrePunto), ''; ...
    'Modelo ajustado (menor SSR ponderado)', nombreModelo(m.tipo), ''; ...
    'Pepita C0 (nugget)', sprintf('%.6g', m.nugget), 'm^2'; ...
    'Meseta parcial C1 (sill)', sprintf('%.6g', m.sill), 'm^2'; ...
    'Meseta total C0 + C1 = C(0)', sprintf('%.6g', cTot), 'm^2'; ...
    'Rango a (parametro del modelo)', sprintf('%.4f', m.rango), 'km'; ...
    'Rango practico (h con 95 % de la meseta parcial)', sprintf('%.4f', rPrac), 'km'; ...
    'SSR ponderado del ajuste (MCP, peso = no. de pares)', sprintf('%.6g', m.SSRw), '-'; ...
    'Distancia maxima de rezago (lagMax = max(d)/2)', sprintf('%.4f', m.lagMax), 'km'; ...
    'Rezagos con pares / configurados', sprintf('%d / %d', numel(m.hExp), par.nLags), '-'; ...
    'Total de pares usados en el variograma experimental', sprintf('%d', sum(m.nPares)), '-'; ...
    'Nodos en la ventana de ajuste', sprintf('%d', m.nNodos), '-'; ...
    'Semiancho efectivo de la ventana', sprintf('%d', semiEfectivo), 'nodos'; ...
    'Vecinos en el sistema de kriging', sprintf('%d', nVec), '-'; ...
    'sigma_KO (desv. estandar de kriging)', sprintf('%.4f', r.sigma), 'm'; ...
    'Semivariograma / covariograma', 'gamma(0) = 0; C(h) = (C0 + C1) - gamma(h) para h > 0', ''; ...
    'Alcance del ajuste', 'Variograma local: se reajusta para cada punto consultado', ''};
end
% [FIN NUEVO]