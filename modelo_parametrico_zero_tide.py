#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modelo paramétrico GNSS/nivelación en Zero Tide con 6 métodos de interpolación.

Flujo
-----
1. Recepción y verificación del sistema de mareas.
   - Se lee el archivo de puntos (lat, lon, h y/o H). El usuario indica qué
     columna es cada cosa (por nombre o número) y en qué sistema de mareas
     viene cada altura y cada malla de N (zero_tide, mean_tide o tide_free).
   - Fórmulas del paper (Guimarães et al. 2025, Earth Planets Space 77:98, §2.4):
       Ec. (1)  H_ZT = H_MT + ΔW_ZT/g
       Ec. (2)  ΔW_ZT/g ≈ -0.198 (3/2 sin²φ - 1/2)
       Ec. (3)  h_ZT = h_TF + 0.0602 - 0.1790 sin²φ - 0.0018 sin⁴φ
     Lo que ya viene en zero tide NO se transforma.
   - Depuración de puntos con errores (reporte en puntos_rechazados.csv).
2. Exportación de las coordenadas en Zero Tide a CSV.
3. Interpolación de la(s) malla(s) de N (geoide/cuasigeoide) en los puntos
   con los 6 métodos de appQGeoidCOL2023_0.m: vecino más cercano, lineal por
   triángulos, bilineal, bicúbica (Keys a=-0.5), IDW y kriging ordinario
   (variograma local ajustado por mínimos cuadrados ponderados).
4. Modelo paramétrico del paper, Ec. (18), resuelto por mínimos cuadrados con
   P = I (Ec. 19-20):
       l_P = h_P - H_P - N_P = δH + x1 (φ_P - φ0) + x2 (λ_P - λ0) cos φ_P
   Para cada método se calcula R² ajustado; el de mayor R²_adj es el mejor.
   (Opcional: --modelo 23 usa la superficie de corrección de 4 parámetros.)

Solo usa la librería estándar de Python (>= 3.8); no requiere numpy.

Uso
---
  python modelo_parametrico_zero_tide.py                  # interactivo
  python modelo_parametrico_zero_tide.py --no-interactivo \\
      --marea-h zero_tide --marea-H zero_tide \\
      --marea-malla zero_tide --marea-malla zero_tide
Ejecute con --help para ver todas las opciones.
"""

import argparse
import bisect
import csv
import math
import os
import statistics
import sys
import time

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------
SISTEMAS = ('zero_tide', 'mean_tide', 'tide_free')
ALIAS_SISTEMA = {
    'zt': 'zero_tide', 'zero': 'zero_tide', 'zero_tide': 'zero_tide', 'zero-tide': 'zero_tide',
    'mt': 'mean_tide', 'mean': 'mean_tide', 'mean_tide': 'mean_tide', 'mean-tide': 'mean_tide',
    'tf': 'tide_free', 'free': 'tide_free', 'tide_free': 'tide_free', 'tide-free': 'tide_free',
    'non_tidal': 'tide_free', 'non-tidal': 'tide_free',
}
K_LOVE = 0.3          # número de Love k (Ekman 1989), solo para N y H en tide free
EPS = 2.220446049250313e-16
METODOS = ('Vecino mas cercano', 'Lineal por triangulos', 'Bilineal', 'Bicubica',
           'IDW', 'Kriging ordinario')
CLAVE_METODO = {'Vecino mas cercano': 'vecino', 'Lineal por triangulos': 'triangulos',
                'Bilineal': 'bilineal', 'Bicubica': 'bicubica', 'IDW': 'idw',
                'Kriging ordinario': 'kriging'}
BITACORA = []


def log(texto=''):
    print(texto)
    BITACORA.append(texto)


# ---------------------------------------------------------------------------
# 1. Sistema de mareas (paper §2.4)
# ---------------------------------------------------------------------------
def delta_H_mt_a_zt(lat):
    """Ec. (1)-(2): H_ZT - H_MT = ΔW_ZT/g ≈ -0.198 (3/2 sin²φ - 1/2) [m]."""
    s2 = math.sin(math.radians(lat)) ** 2
    return -0.198 * (1.5 * s2 - 0.5)


def delta_h_tf_a_zt(lat):
    """Ec. (3): h_ZT - h_TF = 0.0602 - 0.1790 sin²φ - 0.0018 sin⁴φ [m]."""
    s2 = math.sin(math.radians(lat)) ** 2
    return 0.0602 - 0.1790 * s2 - 0.0018 * s2 * s2


def correccion_h(lat, sistema):
    """Corrección a sumar a h para llevarla a zero tide."""
    if sistema == 'zero_tide':
        return 0.0
    if sistema == 'tide_free':
        return delta_h_tf_a_zt(lat)                       # Ec. (3)
    # mean_tide: la corteza 'mean tide' coincide con la 'zero tide' (la
    # deformación permanente está incluida en ambas, Ekman 1989).
    return 0.0


def correccion_N(lat, sistema):
    """Corrección a sumar a N (o ζ) para llevarla a zero tide.

    No está escrita en el paper; se deriva de forma consistente con la Ec. (2):
    N_MT - N_ZT = δ(φ) = -0.198(3/2 sin²φ - 1/2) y N_ZT - N_TF = k·δ(φ)
    (Ekman 1989; Mäkinen & Ihde 2009), con k = 0.3.
    """
    if sistema == 'zero_tide':
        return 0.0
    d = delta_H_mt_a_zt(lat)
    if sistema == 'mean_tide':
        return -d
    return K_LOVE * d


def correccion_H(lat, sistema):
    """Corrección a sumar a H para llevarla a zero tide."""
    if sistema == 'zero_tide':
        return 0.0
    if sistema == 'mean_tide':
        return delta_H_mt_a_zt(lat)                       # Ec. (1)-(2)
    # tide_free (poco habitual): H = h - N -> ΔH = Δh(Ec. 3) - ΔN
    return delta_h_tf_a_zt(lat) - correccion_N(lat, 'tide_free')


def normalizar_sistema(txt):
    if txt is None:
        return None
    return ALIAS_SISTEMA.get(str(txt).strip().lower())


# ---------------------------------------------------------------------------
# Interacción con el usuario
# ---------------------------------------------------------------------------
def es_interactivo(args):
    return not args.no_interactivo and sys.stdin.isatty()


def preguntar(args, mensaje, defecto=None):
    if not es_interactivo(args):
        return defecto
    sufijo = f' [{defecto}]' if defecto not in (None, '') else ''
    r = input(f'{mensaje}{sufijo}: ').strip()
    return r if r else defecto


def preguntar_sistema(args, que, dado, sugerido=None):
    """Devuelve el sistema de mareas de 'que' (pregunta si no se dio)."""
    sis = normalizar_sistema(dado)
    if dado is not None and sis is None:
        sys.exit(f'ERROR: sistema de mareas no reconocido para {que}: "{dado}". '
                 f'Use uno de {SISTEMAS}.')
    if sis is None and es_interactivo(args):
        while sis is None:
            r = preguntar(args, f'¿En qué sistema de mareas está {que}? '
                                f'(zero_tide / mean_tide / tide_free)', sugerido)
            sis = normalizar_sistema(r)
            if sis is None:
                print('  Respuesta no válida. Escriba zero_tide, mean_tide o tide_free (o zt/mt/tf).')
    if sis is None:
        sys.exit(f'ERROR: no se indicó el sistema de mareas de {que}. En modo no interactivo '
                 f'use las opciones --marea-h / --marea-H / --marea-malla.')
    return sis


# ---------------------------------------------------------------------------
# Lectura y validación de mallas (adaptado de cargarYValidarMalla del .m)
# ---------------------------------------------------------------------------
class Malla:
    """Malla regular: lat y lon ascendentes, Z[i][j] (None si falta el nodo)."""

    def __init__(self, nombre, lat, lon, Z, sistema_cabecera=None):
        self.nombre = nombre
        self.lat = lat
        self.lon = lon
        self.Z = Z
        self.sistema_cabecera = sistema_cabecera
        self.sistema = None
        self.xy = None          # proyección Origen-Nacional de los nodos (km)

    def preparar_proyeccion(self):
        self.xy = [[proyectar_origen_nacional(la, lo) for lo in self.lon] for la in self.lat]


def _valores_unicos(v, tol):
    v = sorted(v)
    u = [v[0]]
    for x in v[1:]:
        if x - u[-1] > tol:
            u.append(x)
    return u


def cargar_malla(ruta, orden=None):
    """Lee una malla 'lat lon N' (sin encabezado) o un archivo ICGEM (.gdf)."""
    if not os.path.isfile(ruta):
        sys.exit(f'ERROR: no se encontró la malla "{ruta}".')
    nombre = os.path.splitext(os.path.basename(ruta))[0]
    with open(ruta, encoding='utf-8', errors='replace') as f:
        lineas = f.read().splitlines()

    cabecera = {}
    inicio = 0
    for k, l in enumerate(lineas[:200]):
        if l.strip().startswith('end_of_head'):
            inicio = k + 1
            for c in lineas[:k]:
                partes = c.split()
                if len(partes) >= 2:
                    cabecera[partes[0].lower()] = partes[1]
            break
    es_icgem = inicio > 0
    if orden is None:
        orden = 'lon,lat' if es_icgem else 'lat,lon'

    datos, n_no_num, n_no_fin = [], 0, 0
    for l in lineas[inicio:]:
        p = l.replace(',', ' ').split()
        if not p:
            continue
        try:
            v = [float(x) for x in p]
        except ValueError:
            n_no_num += 1
            continue
        if len(v) != 3:
            n_no_num += 1
            continue
        if not all(math.isfinite(x) for x in v):
            n_no_fin += 1
            continue
        datos.append((v[1], v[0], v[2]) if orden == 'lon,lat' else tuple(v))

    gap = float(cabecera['gapvalue']) if 'gapvalue' in cabecera else None
    n_gap = 0
    if gap is not None:
        n_gap = sum(1 for d in datos if abs(d[2] - gap) < 1e-9)
        datos = [d for d in datos if abs(d[2] - gap) >= 1e-9]
    # Longitudes 0..360 -> -180..180
    datos = [(la, lo - 360.0 if lo > 180.0 else lo, z) for la, lo, z in datos]
    fuera = [d for d in datos if abs(d[0]) > 90 or abs(d[1]) > 180]
    datos = [d for d in datos if abs(d[0]) <= 90 and abs(d[1]) <= 180]
    if len(datos) < 16:
        sys.exit(f'ERROR: la malla "{ruta}" tiene {len(datos)} registros válidos.')

    latU = _valores_unicos([d[0] for d in datos], 1e-6)
    lonU = _valores_unicos([d[1] for d in datos], 1e-6)
    nLat, nLon = len(latU), len(lonU)
    dLat = (latU[-1] - latU[0]) / (nLat - 1)
    dLon = (lonU[-1] - lonU[0]) / (nLon - 1)
    desvLat = max(abs(b - a - dLat) for a, b in zip(latU, latU[1:]))
    desvLon = max(abs(b - a - dLon) for a, b in zip(lonU, lonU[1:]))
    if desvLat > 0.01 * dLat or desvLon > 0.01 * dLon:
        sys.exit(f'ERROR: la malla "{ruta}" no es regular (desviación del paso '
                 f'{desvLat:.3g}° en lat, {desvLon:.3g}° en lon).')
    lat = [latU[0] + k * dLat for k in range(nLat)]
    lon = [lonU[0] + k * dLon for k in range(nLon)]
    Z = [[None] * nLon for _ in range(nLat)]
    n_dup = n_conf = 0
    for la, lo, z in datos:
        i = round((la - lat[0]) / dLat)
        j = round((lo - lon[0]) / dLon)
        if Z[i][j] is not None:
            n_dup += 1
            if abs(Z[i][j] - z) > 1e-6:
                n_conf += 1
            continue
        Z[i][j] = z
    n_falt = sum(1 for fila in Z for z in fila if z is None)
    vals = [z for fila in Z for z in fila if z is not None]

    m = Malla(nombre, lat, lon, Z, normalizar_sistema(cabecera.get('tide_system')))
    log(f'  Malla {nombre} ({"ICGEM" if es_icgem else "lat lon N"}; columnas {orden})')
    log(f'    nodos válidos {len(datos)}, no numéricos {n_no_num}, no finitos {n_no_fin}, '
        f'gapvalue {n_gap}, fuera de rango {len(fuera)}')
    log(f'    {nLat} x {nLon} nodos, paso {dLat * 60:.4f}\' x {dLon * 60:.4f}\', '
        f'lat [{lat[0]:.4f}, {lat[-1]:.4f}], lon [{lon[0]:.4f}, {lon[-1]:.4f}]')
    log(f'    duplicados {n_dup} (con valor distinto {n_conf}), nodos faltantes {n_falt}, '
        f'N en [{min(vals):.3f}, {max(vals):.3f}] m')
    if cabecera.get('tide_system'):
        log(f'    El encabezado declara tide_system = {cabecera["tide_system"]}')
    if any(abs(z) > 150 for z in vals):
        log('    ADVERTENCIA: hay valores |N| > 150 m (no plausibles).')
    return m


def aplicar_marea_malla(m):
    """Lleva los nodos de la malla a zero tide (si hace falta)."""
    if m.sistema == 'zero_tide':
        log(f'  Malla {m.nombre}: zero_tide -> sin transformación.')
        return
    for i, la in enumerate(m.lat):
        c = correccion_N(la, m.sistema)
        fila = m.Z[i]
        for j in range(len(fila)):
            if fila[j] is not None:
                fila[j] += c
    log(f'  Malla {m.nombre}: {m.sistema} -> zero_tide '
        f'(N_ZT = N + corrección; derivada de la Ec. 2, k = {K_LOVE}).')


# ---------------------------------------------------------------------------
# Proyección MAGNA-SIRGAS / Origen-Nacional (EPSG:9377), igual que en el .m
# ---------------------------------------------------------------------------
_A = 6378137.0
_F = 1 / 298.257222101
_E2 = 2 * _F - _F ** 2
_EP2 = _E2 / (1 - _E2)


def _arco_meridiano(p):
    e2 = _E2
    return _A * ((1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256) * p
                 - (3 * e2 / 8 + 3 * e2 ** 2 / 32 + 45 * e2 ** 3 / 1024) * math.sin(2 * p)
                 + (15 * e2 ** 2 / 256 + 45 * e2 ** 3 / 1024) * math.sin(4 * p)
                 - (35 * e2 ** 3 / 3072) * math.sin(6 * p))


_M0 = _arco_meridiano(math.radians(4.0))


def proyectar_origen_nacional(lat, lon):
    """Devuelve (x, y) en km."""
    k0, FE, FN = 0.9992, 5000000.0, 2000000.0
    phi = math.radians(lat)
    lam = math.radians(lon)
    N = _A / math.sqrt(1 - _E2 * math.sin(phi) ** 2)
    T = math.tan(phi) ** 2
    C = _EP2 * math.cos(phi) ** 2
    A = (lam - math.radians(-73.0)) * math.cos(phi)
    x = FE + k0 * N * (A + (1 - T + C) * A ** 3 / 6
                       + (5 - 18 * T + T ** 2 + 72 * C - 58 * _EP2) * A ** 5 / 120)
    y = FN + k0 * (_arco_meridiano(phi) - _M0 + N * math.tan(phi) * (
        A ** 2 / 2 + (5 - T + 9 * C + 4 * C ** 2) * A ** 4 / 24
        + (61 - 58 * T + T ** 2 + 600 * C - 330 * _EP2) * A ** 6 / 720))
    return x / 1000.0, y / 1000.0


# ---------------------------------------------------------------------------
# 3. Seis métodos de interpolación (emulan appQGeoidCOL2023_0.m)
#    Cada método devuelve (nodos [(i, j)], pesos) y N = Σ w_i N_i.
# ---------------------------------------------------------------------------
class NoAplicable(Exception):
    pass


def _celda_1d(v, x):
    tol = 1e-9
    if x < v[0] - tol or x > v[-1] + tol:
        return None
    i = bisect.bisect_right(v, x + tol) - 1
    return min(max(i, 0), len(v) - 2)


def ubicar_celda(m, lat0, lon0):
    i = _celda_1d(m.lat, lat0)
    j = _celda_1d(m.lon, lon0)
    if i is None or j is None:
        raise NoAplicable('fuera de la cobertura de la malla (no se extrapola)')
    v = min(max((lat0 - m.lat[i]) / (m.lat[i + 1] - m.lat[i]), 0.0), 1.0)
    u = min(max((lon0 - m.lon[j]) / (m.lon[j + 1] - m.lon[j]), 0.0), 1.0)
    return i, j, u, v


def ventana(m, i, j, r):
    """Nodos (no faltantes) de una ventana 2r x 2r alrededor de la celda (i, j),
    en el mismo orden que meshgrid(...)(:) de MATLAB (columna a columna)."""
    ii = range(max(0, i - r + 1), min(len(m.lat) - 1, i + r) + 1)
    jj = range(max(0, j - r + 1), min(len(m.lon) - 1, j + r) + 1)
    return [(a, b) for b in jj for a in ii if m.Z[a][b] is not None]


def pesos_vecino(m, lat0, lon0, xy0):
    i, j, _, _ = ubicar_celda(m, lat0, lon0)
    esquinas = [(i, j), (i, j + 1), (i + 1, j + 1), (i + 1, j)]   # SW, SE, NE, NW
    d = [math.dist(m.xy[a][b], xy0) for a, b in esquinas]
    k = d.index(min(d))
    return [esquinas[k]], [1.0]


def pesos_triangulos(m, lat0, lon0, xy0):
    i, j, u, v = ubicar_celda(m, lat0, lon0)
    if u >= v:      # triángulo SW-SE-NE (diagonal SW-NE)
        return [(i, j), (i, j + 1), (i + 1, j + 1)], [1 - u, u - v, v]
    return [(i, j), (i + 1, j + 1), (i + 1, j)], [1 - v, u, v - u]


def pesos_bilineal(m, lat0, lon0, xy0):
    i, j, u, v = ubicar_celda(m, lat0, lon0)
    return ([(i, j), (i, j + 1), (i + 1, j + 1), (i + 1, j)],
            [(1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v])


def _nucleo_keys(s, a=-0.5):
    s = abs(s)
    if s <= 1:
        return (a + 2) * s ** 3 - (a + 3) * s ** 2 + 1
    if s < 2:
        return a * s ** 3 - 5 * a * s ** 2 + 8 * a * s - 4 * a
    return 0.0


def pesos_bicubica(m, lat0, lon0, xy0):
    i, j, u, v = ubicar_celda(m, lat0, lon0)
    if i - 1 < 0 or i + 2 > len(m.lat) - 1 or j - 1 < 0 or j + 2 > len(m.lon) - 1:
        raise NoAplicable('soporte 4x4 incompleto (celda de borde)')
    wx = [_nucleo_keys(s) for s in (u + 1, u, 1 - u, 2 - u)]
    wy = [_nucleo_keys(s) for s in (v + 1, v, 1 - v, 2 - v)]
    nodos, w = [], []
    for b in range(4):
        for a in range(4):
            nodos.append((i - 1 + a, j - 1 + b))
            w.append(wy[a] * wx[b])
    return nodos, w


def pesos_idw(m, lat0, lon0, xy0, p=2.0, k=16):
    i, j, _, _ = ubicar_celda(m, lat0, lon0)
    r = math.ceil(math.sqrt(k) / 2) + 1
    cand = ventana(m, i, j, r)
    d = [math.dist(m.xy[a][b], xy0) for a, b in cand]
    orden = sorted(range(len(d)), key=d.__getitem__)[:min(k, len(d))]
    nodos = [cand[o] for o in orden]
    dk = [d[o] for o in orden]
    if dk[0] < 1e-6:          # coincide con un nodo (< 1 mm)
        return nodos, [1.0] + [0.0] * (len(dk) - 1)
    w = [1.0 / x ** p for x in dk]
    s = sum(w)
    return nodos, [x / s for x in w]


# --- Kriging ordinario con variograma local ajustado por MCP --------------
def modelo_variograma(h, nugget, sill, rango, tipo):
    if h <= 1e-12:
        return 0.0
    if tipo == 'sph':
        r = min(h / rango, 1.0)
        return nugget + sill * (1.5 * r - 0.5 * r ** 3)
    if tipo == 'exp':
        return nugget + sill * (1 - math.exp(-h / rango))
    if tipo == 'gau':
        return nugget + sill * (1 - math.exp(-(h / rango) ** 2))
    if tipo == 'cir':
        r = min(h / rango, 1.0)
        return nugget + sill * (1 - (2 / math.pi) * (math.acos(r) - r * math.sqrt(1 - r * r)))
    raise ValueError(tipo)


def _eps_de(x):
    return math.ulp(abs(x)) if x else EPS


def fminsearch(f, x0, tolx=1e-4, tolf=1e-4):
    """Nelder-Mead con los mismos parámetros y criterios que fminsearch de MATLAB."""
    n = len(x0)
    maxit = maxfun = 200 * n
    rho, chi, psi, sig = 1.0, 2.0, 0.5, 0.5
    v = [list(x0)]
    fv = [f(x0)]
    for k in range(n):
        y = list(x0)
        y[k] = (1 + 0.05) * y[k] if y[k] != 0 else 0.00025
        v.append(y)
        fv.append(f(y))
    orden = sorted(range(n + 1), key=fv.__getitem__)
    v = [v[o] for o in orden]
    fv = [fv[o] for o in orden]
    it, nf = 1, n + 1
    while nf < maxfun and it < maxit:
        if (max(abs(fv[0] - fv[k]) for k in range(1, n + 1)) <= max(tolf, 10 * _eps_de(fv[0])) and
                max(abs(v[k][c] - v[0][c]) for k in range(1, n + 1) for c in range(n))
                <= max(tolx, 10 * _eps_de(max(abs(c) for c in v[0])))):
            break
        xbar = [sum(v[k][c] for k in range(n)) / n for c in range(n)]
        xr = [(1 + rho) * xbar[c] - rho * v[n][c] for c in range(n)]
        fxr = f(xr)
        nf += 1
        encoger = False
        if fxr < fv[0]:
            xe = [(1 + rho * chi) * xbar[c] - rho * chi * v[n][c] for c in range(n)]
            fxe = f(xe)
            nf += 1
            v[n], fv[n] = (xe, fxe) if fxe < fxr else (xr, fxr)
        elif fxr < fv[n - 1]:
            v[n], fv[n] = xr, fxr
        else:
            if fxr < fv[n]:
                xc = [(1 + psi * rho) * xbar[c] - psi * rho * v[n][c] for c in range(n)]
                fxc = f(xc)
                nf += 1
                if fxc <= fxr:
                    v[n], fv[n] = xc, fxc
                else:
                    encoger = True
            else:
                xcc = [(1 - psi) * xbar[c] + psi * v[n][c] for c in range(n)]
                fxcc = f(xcc)
                nf += 1
                if fxcc < fv[n]:
                    v[n], fv[n] = xcc, fxcc
                else:
                    encoger = True
            if encoger:
                for k in range(1, n + 1):
                    v[k] = [v[0][c] + sig * (v[k][c] - v[0][c]) for c in range(n)]
                    fv[k] = f(v[k])
                nf += n
        orden = sorted(range(n + 1), key=fv.__getitem__)
        v = [v[o] for o in orden]
        fv = [fv[o] for o in orden]
        it += 1
    return v[0]


def ajustar_variograma(xy, z, n_lags=12, tipos=('sph', 'exp', 'gau', 'cir')):
    """Semivariograma experimental de la ventana + ajuste MCP (peso = nº de pares)."""
    n = len(z)
    if n < 6:
        raise NoAplicable(f'solo {n} nodos en la ventana para el variograma')
    dist, semiv = [], []
    for a in range(n):
        xa, ya = xy[a]
        za = z[a]
        for b in range(a + 1, n):
            dist.append(math.hypot(xa - xy[b][0], ya - xy[b][1]))
            semiv.append(0.5 * (za - z[b]) ** 2)
    lag_max = max(dist) / 2
    bordes = [lag_max * k / n_lags for k in range(n_lags + 1)]
    sh = [0.0] * n_lags
    sg = [0.0] * n_lags
    npar = [0] * n_lags
    for d, s in zip(dist, semiv):
        k = bisect.bisect_left(bordes, d) - 1      # bordes[k] < d <= bordes[k+1]
        if 0 <= k < n_lags:
            sh[k] += d
            sg[k] += s
            npar[k] += 1
    h = [sh[k] / npar[k] for k in range(n_lags) if npar[k]]
    g = [sg[k] / npar[k] for k in range(n_lags) if npar[k]]
    w = [npar[k] for k in range(n_lags) if npar[k]]
    if len(h) < 3:
        raise NoAplicable('menos de 3 rezagos con pares')
    ini = [math.sqrt(max(min(g) * 0.3, EPS)), math.sqrt(max(max(g) - min(g), EPS)), 0.0]

    def rango(p3):
        return lag_max / (1 + math.exp(-max(min(p3, 700), -700))) + EPS

    mejor = None
    for tipo in tipos:
        def obj(p, tipo=tipo):
            nug, sil, rg = p[0] ** 2, p[1] ** 2, rango(p[2])
            return sum(wk * (gk - modelo_variograma(hk, nug, sil, rg, tipo)) ** 2
                       for hk, gk, wk in zip(h, g, w))
        p = fminsearch(obj, ini)
        ssr = obj(p)
        if mejor is None or ssr < mejor['SSRw']:
            mejor = {'tipo': tipo, 'nugget': p[0] ** 2, 'sill': p[1] ** 2,
                     'rango': rango(p[2]), 'SSRw': ssr}
    return mejor


def resolver(A, b):
    """Eliminación gaussiana con pivoteo parcial. Devuelve (x, pivote_min_relativo)."""
    n = len(b)
    M = [fila[:] + [b[k]] for k, fila in enumerate(A)]
    escala = max(abs(x) for fila in A for x in fila) or 1.0
    piv_min = float('inf')
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-300:
            raise NoAplicable('sistema singular')
        M[c], M[p] = M[p], M[c]
        piv_min = min(piv_min, abs(M[c][c]) / escala)
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f:
                for k in range(c, n + 1):
                    M[r][k] -= f * M[c][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (M[r][n] - sum(M[r][k] * x[k] for k in range(r + 1, n))) / M[r][r]
    return x, piv_min


def pesos_kriging(m, lat0, lon0, xy0, n_krig=16, semiancho=6, n_lags=12, cache=None):
    i, j, _, _ = ubicar_celda(m, lat0, lon0)
    r = max(semiancho, math.ceil(math.sqrt(n_krig) / 2) + 1)
    cand = ventana(m, i, j, r)
    xy = [m.xy[a][b] for a, b in cand]
    clave = (i, j)
    if cache is not None and clave in cache:
        mod = cache[clave]
    else:
        mod = ajustar_variograma(xy, [m.Z[a][b] for a, b in cand], n_lags)
        if cache is not None:
            cache[clave] = mod
    d = [math.dist(p, xy0) for p in xy]
    orden = sorted(range(len(d)), key=d.__getitem__)[:min(n_krig, len(d))]
    k = len(orden)
    gm = lambda hh: modelo_variograma(hh, mod['nugget'], mod['sill'], mod['rango'], mod['tipo'])
    A = [[gm(math.dist(xy[orden[a]], xy[orden[b]])) for b in range(k)] + [1.0] for a in range(k)]
    A.append([1.0] * k + [0.0])
    g0 = [gm(d[o]) for o in orden]
    sol, piv = resolver(A, g0 + [1.0])
    lam, mu = sol[:k], sol[k]
    var = sum(l * g for l, g in zip(lam, g0)) + mu
    info = {'sigma': math.sqrt(max(var, 0.0)), 'modelo': mod['tipo'], 'mal_cond': piv < 1e-14}
    return [cand[o] for o in orden], lam, info


def interpolar(m, metodo, lat0, lon0, xy0, par, cache):
    """Devuelve (N, suma_pesos, info) o lanza NoAplicable."""
    info = {}
    if metodo == 'Vecino mas cercano':
        nodos, w = pesos_vecino(m, lat0, lon0, xy0)
    elif metodo == 'Lineal por triangulos':
        nodos, w = pesos_triangulos(m, lat0, lon0, xy0)
    elif metodo == 'Bilineal':
        nodos, w = pesos_bilineal(m, lat0, lon0, xy0)
    elif metodo == 'Bicubica':
        nodos, w = pesos_bicubica(m, lat0, lon0, xy0)
    elif metodo == 'IDW':
        nodos, w = pesos_idw(m, lat0, lon0, xy0, par.idw_p, par.idw_k)
    else:
        nodos, w, info = pesos_kriging(m, lat0, lon0, xy0, par.krig_n, par.krig_ventana,
                                       12, cache)
    z = 0.0
    for (a, b), wi in zip(nodos, w):
        zi = m.Z[a][b]
        if zi is None:
            if wi != 0:
                raise NoAplicable('falta un nodo de la malla con peso no nulo')
            continue
        z += wi * zi
    return z, sum(w), info


# ---------------------------------------------------------------------------
# 4. Modelo paramétrico (Ec. 18-20) y R² ajustado
# ---------------------------------------------------------------------------
def invertir(M):
    n = len(M)
    inv = []
    for c in range(n):
        e = [1.0 if r == c else 0.0 for r in range(n)]
        x, _ = resolver(M, e)
        inv.append(x)
    return [[inv[c][r] for c in range(n)] for r in range(n)]


def ajustar_modelo(puntos, l, modelo='18'):
    """Mínimos cuadrados x = (AᵀPA)⁻¹AᵀPl con P = I."""
    n = len(l)
    phi0 = statistics.fmean(p['lat'] for p in puntos)
    lam0 = statistics.fmean(p['lon'] for p in puntos)
    filas = []
    for p in puntos:
        f, la = p['lat'], p['lon']
        if modelo == '18':     # δH + x1(φ-φ0) + x2(λ-λ0)cosφ   [φ, λ en grados]
            filas.append([1.0, f - phi0, (la - lam0) * math.cos(math.radians(f))])
        else:                  # x0 + x1 cosφcosλ + x2 cosφsinλ + x3 sinφ  (Ec. 23)
            fr, lr = math.radians(f), math.radians(la)
            filas.append([1.0, math.cos(fr) * math.cos(lr), math.cos(fr) * math.sin(lr),
                          math.sin(fr)])
    u = len(filas[0])
    N = [[sum(r[a] * r[b] for r in filas) for b in range(u)] for a in range(u)]
    t = [sum(r[a] * li for r, li in zip(filas, l)) for a in range(u)]
    x, _ = resolver(N, t)
    v = [li - sum(ra * xa for ra, xa in zip(r, x)) for r, li in zip(filas, l)]
    sse = sum(e * e for e in v)
    media = statistics.fmean(l)
    sst = sum((li - media) ** 2 for li in l)
    r2 = 1 - sse / sst if sst > 0 else float('nan')
    r2adj = 1 - (1 - r2) * (n - 1) / (n - u)
    s02 = sse / (n - u)
    Qx = invertir(N)
    sx = [math.sqrt(s02 * Qx[k][k]) for k in range(u)]
    return {'x': x, 'sx': sx, 'v': v, 'n': n, 'u': u, 'R2': r2, 'R2adj': r2adj,
            'sigma0': math.sqrt(s02), 'rmse': math.sqrt(sse / n),
            'media_l': media, 'std_l': statistics.stdev(l), 'phi0': phi0, 'lam0': lam0}


# ---------------------------------------------------------------------------
# Lectura y depuración de puntos
# ---------------------------------------------------------------------------
def leer_tabla(ruta):
    if not os.path.isfile(ruta):
        sys.exit(f'ERROR: no se encontró el archivo de puntos "{ruta}".')
    for enc in ('utf-8-sig', 'cp1252', 'latin-1'):
        try:
            with open(ruta, encoding=enc) as f:
                texto = f.read()
            break
        except UnicodeDecodeError:
            continue
    lineas = [l for l in texto.splitlines() if l.strip()]
    delim = None
    for d in (';', '\t', ','):
        if all(l.count(d) == lineas[0].count(d) and l.count(d) > 0 for l in lineas[:20]):
            delim = d
            break
    if delim is None:
        filas = [l.split() for l in lineas]
    else:
        filas = list(csv.reader(lineas, delimiter=delim))
    try:                       # ¿la primera fila es numérica? -> sin encabezado
        [float(x) for x in filas[0]]
        encabezado = [f'col{k + 1}' for k in range(len(filas[0]))]
    except ValueError:
        encabezado = [c.strip() for c in filas[0]]
        filas = filas[1:]
    return encabezado, filas, delim


SUGERENCIAS = {
    'id': ('fid', 'id', 'punto', 'nombre', 'name', 'codigo'),
    'lat': ('latitud', 'lat', 'latitude', 'phi'),
    'lon': ('longitud', 'lon', 'long', 'longitude', 'lambda'),
    'h': ('alt_elipso', 'h', 'h_elip', 'altura_elipsoidal', 'ellipsoidal_height', 'elip'),
    'H': ('altura_m_s', 'h_nivel', 'altura_nivelada', 'ortometrica', 'altura_msnm', 'msnm',
          'orthometric_height', 'h_orto', 'hn'),
}


def sugerir_columna(encabezado, clave):
    low = [c.lower() for c in encabezado]
    for s in SUGERENCIAS[clave]:
        if s in low and not (clave == 'h' and encabezado[low.index(s)] == 'H'):
            return encabezado[low.index(s)]
    return None


def resolver_columna(encabezado, valor):
    if valor in (None, '', '-', 'ninguna', 'no'):
        return None
    if valor in encabezado:
        return encabezado.index(valor)
    low = [c.lower() for c in encabezado]
    if str(valor).lower() in low:
        return low.index(str(valor).lower())
    try:
        k = int(valor) - 1
        if 0 <= k < len(encabezado):
            return k
    except ValueError:
        pass
    sys.exit(f'ERROR: la columna "{valor}" no existe. Columnas: {encabezado}')


def mapear_columnas(args, encabezado):
    log('\nColumnas del archivo de puntos:')
    for k, c in enumerate(encabezado, 1):
        log(f'  {k}: {c}')
    textos = {'id': 'identificador del punto (opcional)', 'lat': 'LATITUD (grados)',
              'lon': 'LONGITUD (grados)', 'h': 'ALTURA ELIPSOIDAL h (m) — vacío/"-" si no viene',
              'H': 'ALTURA SOBRE EL NIVEL DEL MAR H (nivelada, m) — vacío/"-" si no viene'}
    dados = {'id': args.col_id, 'lat': args.col_lat, 'lon': args.col_lon, 'h': args.col_h,
             'H': args.col_H}
    idx = {}
    for clave in ('id', 'lat', 'lon', 'h', 'H'):
        val = dados[clave]
        if val is None:
            val = preguntar(args, f'Columna de {textos[clave]}', sugerir_columna(encabezado, clave))
        idx[clave] = resolver_columna(encabezado, val)
    if idx['lat'] is None or idx['lon'] is None:
        sys.exit('ERROR: hay que indicar las columnas de latitud y longitud.')
    if idx['h'] is None and idx['H'] is None:
        sys.exit('ERROR: hay que indicar al menos una altura (h elipsoidal o H sobre el nivel del mar).')
    for clave in ('id', 'lat', 'lon', 'h', 'H'):
        log(f'  {clave:>3} -> {encabezado[idx[clave]] if idx[clave] is not None else "(no viene)"}')
    return idx


def _num(txt):
    try:
        v = float(str(txt).strip().replace(',', '.'))
        return v if math.isfinite(v) else None
    except ValueError:
        return None


def leer_puntos(filas, idx, rechazos):
    """Reglas R1-R5 (formato, rango, duplicados, h nula, N=h-H no plausible)."""
    puntos = []
    vistos_id, vistos_xy = {}, {}
    for k, fila in enumerate(filas, start=2):
        pid = fila[idx['id']].strip() if idx['id'] is not None and idx['id'] < len(fila) else str(k - 1)
        get = lambda c: _num(fila[idx[c]]) if idx[c] is not None and idx[c] < len(fila) else None
        lat, lon = get('lat'), get('lon')
        h, H = get('h'), get('H')
        base = {'id': pid, 'linea': k, 'lat': lat, 'lon': lon, 'h': h, 'H': H}

        def rechazar(regla, detalle):
            rechazos.append(dict(base, regla=regla, detalle=detalle))

        faltan = [c for c in ('lat', 'lon', 'h', 'H') if idx[c] is not None and get(c) is None]
        if faltan:
            rechazar('R1 dato vacío/no numérico', f'campos: {", ".join(faltan)}')
            continue
        if lon > 180:
            lon -= 360.0
            base['lon'] = lon
        if abs(lat) > 90 or abs(lon) > 180:
            rechazar('R2 coordenada fuera de rango', f'lat={lat}, lon={lon}')
            continue
        if idx['id'] is not None and pid in vistos_id:
            rechazar('R3 identificador duplicado', f'mismo id que la línea {vistos_id[pid]}')
            continue
        clave_xy = (round(lat, 7), round(lon, 7))
        if clave_xy in vistos_xy:
            rechazar('R3 coordenadas duplicadas', f'mismas lat/lon que la línea {vistos_xy[clave_xy]}')
            continue
        if h is not None and H is not None:
            if abs(h) < 0.01 and abs(H) > 1.0:
                rechazar('R4 altura elipsoidal nula', f'h={h} m con H={H} m')
                continue
            if abs(h - H) > 150.0:
                rechazar('R5 h-H no plausible', f'h-H={h - H:.3f} m (|N| > 150 m)')
                continue
        vistos_id[pid] = k
        vistos_xy[clave_xy] = k
        puntos.append(base)
    return puntos


def depurar_contra_mallas(puntos, mallas, umbral_z, rechazos):
    """R6: error grosero. Con la malla (bilineal) se calcula l = h - H - N; se
    rechazan los puntos con |l - mediana| > umbral_z · 1.4826 · MAD en alguna
    malla. También R7: punto fuera de la cobertura de alguna malla."""
    malos = {}
    for m in mallas:
        ls = []
        for p in puntos:
            try:
                N, _, _ = interpolar(m, 'Bilineal', p['lat'], p['lon'], None, None, None)
            except NoAplicable:
                malos.setdefault(p['id'], ('R7 fuera de cobertura', f'malla {m.nombre}'))
                continue
            ls.append((p, p['h_zt'] - p['H_zt'] - N))
        valores = [x for _, x in ls]
        med = statistics.median(valores)
        mad = statistics.median(abs(x - med) for x in valores)
        lim = umbral_z * 1.4826 * mad
        n_m = 0
        for p, x in ls:
            if abs(x - med) > lim:
                n_m += 1
                if p['id'] not in malos:
                    malos[p['id']] = ('R6 error grosero (h-H-N)',
                                      f'malla {m.nombre}: l={x:.3f} m, mediana={med:.3f} m, '
                                      f'|l-med|={abs(x - med):.3f} > {lim:.3f} m')
        log(f'  R6 con {m.nombre}: mediana(l)={med:.3f} m, σ_robusta={1.4826 * mad:.3f} m, '
            f'límite ±{lim:.3f} m -> {n_m} puntos')
    quedan = []
    for p in puntos:
        if p['id'] in malos:
            regla, det = malos[p['id']]
            rechazos.append({'id': p['id'], 'linea': p['linea'], 'lat': p['lat'], 'lon': p['lon'],
                             'h': p['h'], 'H': p['H'], 'regla': regla, 'detalle': det})
        else:
            quedan.append(p)
    return quedan


# ---------------------------------------------------------------------------
# Utilidades de salida
# ---------------------------------------------------------------------------
def escribir_csv(ruta, campos, filas):
    with open(ruta, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=campos, extrasaction='ignore')
        w.writeheader()
        for r in filas:
            w.writerow({k: ('' if r.get(k) is None else r.get(k)) for k in campos})


def fmt(x, n=4):
    return '' if x is None else f'{x:.{n}f}'


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------
def main():
    aqui = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(
        description='Zero tide + 6 interpolaciones + modelo paramétrico (Ec. 18) con R² ajustado.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--puntos', default=os.path.join(aqui, 'coords_QGEOID2023.txt'),
                    help='archivo de puntos GNSS/nivelación')
    ap.add_argument('--malla', action='append',
                    help='malla de N (repetible). Por defecto QGeoidCOL2023.txt y coords_EGM2008.txt')
    ap.add_argument('--orden-malla', action='append', choices=('lat,lon', 'lon,lat'),
                    help='orden de columnas de cada malla (por defecto: ICGEM=lon,lat; otras=lat,lon)')
    ap.add_argument('--col-id'); ap.add_argument('--col-lat'); ap.add_argument('--col-lon')
    ap.add_argument('--col-h', help='columna de altura elipsoidal ("-" si no viene)')
    ap.add_argument('--col-H', help='columna de altura sobre el nivel del mar ("-" si no viene)')
    ap.add_argument('--marea-h', help='sistema de mareas de h')
    ap.add_argument('--marea-H', help='sistema de mareas de H')
    ap.add_argument('--marea-malla', action='append', help='sistema de mareas de cada malla (repetible)')
    ap.add_argument('--modelo', choices=('18', '23'), default='18',
                    help='18: δH + 2 inclinaciones (paper Ec. 18); 23: superficie de 4 parámetros')
    ap.add_argument('--umbral-mad', type=float, default=5.0,
                    help='R6: rechazo si |l - mediana| > umbral · σ_robusta')
    ap.add_argument('--idw-p', type=float, default=2.0)
    ap.add_argument('--idw-k', type=int, default=16)
    ap.add_argument('--krig-n', type=int, default=16)
    ap.add_argument('--krig-ventana', type=int, default=6, help='semiancho de ventana (nodos)')
    ap.add_argument('--salida', default=os.path.join(aqui, 'resultados'))
    ap.add_argument('--no-interactivo', action='store_true',
                    help='no preguntar: usar opciones de línea de comandos y sugerencias')
    args = ap.parse_args()
    if args.malla is None:
        args.malla = [os.path.join(aqui, f) for f in ('QGeoidCOL2023.txt', 'coords_EGM2008.txt')
                      if os.path.isfile(os.path.join(aqui, f))]
    os.makedirs(args.salida, exist_ok=True)
    t0 = time.time()

    log('=' * 78)
    log('MODELO PARAMÉTRICO GNSS/NIVELACIÓN EN ZERO TIDE (Guimarães et al. 2025)')
    log('=' * 78)

    # ---------------- 1. Recepción, mareas y depuración ----------------
    log(f'\n[1] RECEPCIÓN DE PUNTOS: {args.puntos}')
    encabezado, filas, delim = leer_tabla(args.puntos)
    log(f'  {len(filas)} registros; separador {repr(delim) if delim else "espacios"}')
    idx = mapear_columnas(args, encabezado)
    tiene_h, tiene_H = idx['h'] is not None, idx['H'] is not None

    log('\n  Sistema de mareas de los puntos:')
    sis_h = preguntar_sistema(args, 'la altura elipsoidal h', args.marea_h) if tiene_h else None
    sis_H = preguntar_sistema(args, 'la altura sobre el nivel del mar H', args.marea_H) if tiene_H else None
    if tiene_h:
        log(f'    h: {sis_h} -> ' + ('sin transformación' if sis_h == 'zero_tide' else
            ('Ec. (3)' if sis_h == 'tide_free' else 'sin cambio (corteza mean tide = zero tide)')))
    if tiene_H:
        log(f'    H: {sis_H} -> ' + ('sin transformación' if sis_H == 'zero_tide' else
            ('Ec. (1)-(2)' if sis_H == 'mean_tide' else 'Ec. (3) - ΔN_TF→ZT (derivada)')))

    log('\n[1b] DEPURACIÓN DE PUNTOS (R1-R5)')
    rechazos = []
    puntos = leer_puntos(filas, idx, rechazos)
    for p in puntos:
        p['corr_h'] = correccion_h(p['lat'], sis_h) if tiene_h else None
        p['corr_H'] = correccion_H(p['lat'], sis_H) if tiene_H else None
        p['h_zt'] = p['h'] + p['corr_h'] if tiene_h else None
        p['H_zt'] = p['H'] + p['corr_H'] if tiene_H else None
    log(f'  válidos tras R1-R5: {len(puntos)} de {len(filas)} (rechazados {len(rechazos)})')

    mallas = []
    if tiene_h and tiene_H and args.malla:
        log('\n[1c] MALLAS DE N (geoide / cuasigeoide)')
        for k, ruta in enumerate(args.malla):
            orden = args.orden_malla[k] if args.orden_malla and k < len(args.orden_malla) else None
            m = cargar_malla(ruta, orden)
            dado = args.marea_malla[k] if args.marea_malla and k < len(args.marea_malla) else None
            m.sistema = preguntar_sistema(args, f'la malla {m.nombre}', dado, m.sistema_cabecera)
            if m.sistema_cabecera and m.sistema_cabecera != m.sistema:
                log(f'    AVISO: el encabezado de {m.nombre} dice {m.sistema_cabecera}, pero se '
                    f'indicó {m.sistema}; se usa lo indicado por el usuario.')
            aplicar_marea_malla(m)
            m.preparar_proyeccion()
            mallas.append(m)
        log('\n[1d] DEPURACIÓN CONTRA LAS MALLAS (R6 error grosero, R7 cobertura)')
        puntos = depurar_contra_mallas(puntos, mallas, args.umbral_mad, rechazos)

    conteo = {}
    for r in rechazos:
        conteo[r['regla']] = conteo.get(r['regla'], 0) + 1
    log(f'\n  RESUMEN DE DEPURACIÓN: {len(puntos)} puntos aceptados, {len(rechazos)} rechazados')
    for regla in sorted(conteo):
        log(f'    {regla}: {conteo[regla]}')
    for r in rechazos[:60]:
        log(f'    - id {r["id"]} (línea {r["linea"]}): {r["regla"]} | {r["detalle"]}')
    if len(rechazos) > 60:
        log(f'    ... ({len(rechazos) - 60} más en puntos_rechazados.csv)')
    escribir_csv(os.path.join(args.salida, 'puntos_rechazados.csv'),
                 ['id', 'linea', 'lat', 'lon', 'h', 'H', 'regla', 'detalle'], rechazos)

    # ---------------- 2. Exportación zero tide ----------------
    ruta_zt = os.path.join(args.salida, 'coordenadas_zero_tide.csv')
    filas_zt = [{'id': p['id'], 'lat': fmt(p['lat'], 9), 'lon': fmt(p['lon'], 9),
                 'h_original': fmt(p['h']), 'sistema_h': sis_h, 'corr_h': fmt(p['corr_h'], 5),
                 'h_ZT': fmt(p['h_zt']), 'H_original': fmt(p['H']), 'sistema_H': sis_H,
                 'corr_H': fmt(p['corr_H'], 5), 'H_ZT': fmt(p['H_zt'])} for p in puntos]
    escribir_csv(ruta_zt, ['id', 'lat', 'lon', 'h_original', 'sistema_h', 'corr_h', 'h_ZT',
                           'H_original', 'sistema_H', 'corr_H', 'H_ZT'], filas_zt)
    log(f'\n[2] EXPORTACIÓN: {len(puntos)} puntos en zero tide -> {ruta_zt}')

    if not (tiene_h and tiene_H):
        log('\nEl modelo paramétrico (l = h - H - N) necesita h y H; solo se exportó zero tide.')
        return
    if not mallas:
        log('\nNo se indicó ninguna malla de N; no se puede aplicar el modelo paramétrico.')
        return

    # ---------------- 3 y 4. Interpolación + modelo ----------------
    resumen_global = []
    for m in mallas:
        log('\n' + '=' * 78)
        log(f'[3] INTERPOLACIÓN DE {m.nombre} EN {len(puntos)} PUNTOS (6 métodos)')
        log(f'    IDW p={args.idw_p}, k={args.idw_k}; kriging n={args.krig_n}, '
            f'semiancho={args.krig_ventana} nodos, variograma MCP (sph/exp/gau/cir)')
        cache = {}
        res = {met: {} for met in METODOS}
        no_apl = {met: {} for met in METODOS}
        sum_w = {met: 0.0 for met in METODOS}
        sig_k = {}
        t1 = time.time()
        for n_p, p in enumerate(puntos, 1):
            xy0 = proyectar_origen_nacional(p['lat'], p['lon'])
            for met in METODOS:
                try:
                    z, sw, info = interpolar(m, met, p['lat'], p['lon'], xy0, args, cache)
                    res[met][p['id']] = z
                    sum_w[met] = max(sum_w[met], abs(sw - 1))
                    if info:
                        sig_k[p['id']] = info
                except NoAplicable as e:
                    no_apl[met][p['id']] = str(e)
            if n_p % 500 == 0:
                print(f'    ... {n_p}/{len(puntos)} puntos ({time.time() - t1:.0f} s)', flush=True)
        log(f'    tiempo {time.time() - t1:.1f} s; variogramas ajustados: {len(cache)} celdas')
        for met in METODOS:
            log(f'    {met:<22} aplicable en {len(res[met]):>5} puntos, '
                f'max|Σw-1| = {sum_w[met]:.1e}' +
                (f'; no aplicable: {len(no_apl[met])}' if no_apl[met] else ''))
        mal_cond = sum(1 for i in sig_k.values() if i.get('mal_cond'))
        if mal_cond:
            log(f'    ADVERTENCIA: {mal_cond} sistemas de kriging mal condicionados.')

        comunes = [p for p in puntos if all(p['id'] in res[met] for met in METODOS)]
        if len(comunes) < len(puntos):
            log(f'    Se usan los {len(comunes)} puntos donde los 6 métodos son aplicables '
                f'(mismo conjunto para comparar).')

        log(f'\n[4] MODELO PARAMÉTRICO Ec. ({args.modelo}) sobre l = h_ZT - H_ZT - N_ZT, P = I')
        if args.modelo == '18':
            log('    l = δH + x1(φ-φ0) + x2(λ-λ0)cosφ   (φ, λ en grados; x1, x2 en m/°)')
        else:
            log('    l = x0 + x1 cosφcosλ + x2 cosφsinλ + x3 sinφ')
        ajustes = {}
        for met in METODOS:
            l = [p['h_zt'] - p['H_zt'] - res[met][p['id']] for p in comunes]
            ajustes[met] = ajustar_modelo(comunes, l, args.modelo)
        mejor = max(METODOS, key=lambda mt: ajustes[mt]['R2adj'])

        nombres_x = ['dH', 'x1', 'x2'] if args.modelo == '18' else ['x0', 'x1', 'x2', 'x3']
        cab = f'    {"Método":<22}' + ''.join(f'{c:>16}' for c in nombres_x) + \
              f'{"σ0 [m]":>9}{"R²":>10}{"R²adj":>12}'
        log(cab)
        filas_res = []
        for met in METODOS:
            a = ajustes[met]
            log(f'    {met:<22}' + ''.join(f'{x:>9.4f}±{s:<6.4f}' for x, s in zip(a['x'], a['sx'])) +
                f'{a["sigma0"]:>9.4f}{a["R2"]:>10.6f}{a["R2adj"]:>12.8f}' +
                ('  <-- MEJOR' if met == mejor else ''))
            fila = {'malla': m.nombre, 'metodo': met, 'modelo': f'Ec.{args.modelo}', 'n': a['n'],
                    'parametros': a['u'], 'phi0': fmt(a['phi0'], 6), 'lambda0': fmt(a['lam0'], 6),
                    'media_l': fmt(a['media_l']), 'std_l': fmt(a['std_l']),
                    'sigma0': fmt(a['sigma0']), 'rmse_residuos': fmt(a['rmse']),
                    'R2': fmt(a['R2'], 8), 'R2_ajustado': fmt(a['R2adj'], 8),
                    'mejor': 'SI' if met == mejor else ''}
            for c, x, s in zip(nombres_x, a['x'], a['sx']):
                fila[c] = fmt(x, 6)
                fila['s_' + c] = fmt(s, 6)
            filas_res.append(fila)
            resumen_global.append(fila)
        empates = [mt for mt in METODOS if abs(ajustes[mt]['R2adj'] - ajustes[mejor]['R2adj']) < 1e-10]
        if len(empates) > 1:
            log(f'    Nota: empate en R²adj entre {", ".join(empates)}.')
        log(f'\n    MEJOR MÉTODO para {m.nombre}: {mejor} '
            f'(R²adj = {ajustes[mejor]["R2adj"]:.8f}, σ0 = {ajustes[mejor]["sigma0"]:.4f} m)')
        menor_s0 = min(METODOS, key=lambda mt: ajustes[mt]['sigma0'])
        if menor_s0 != mejor:
            log(f'    (Con el criterio de menor σ0 el mejor sería {menor_s0}.)')

        campos_res = ['malla', 'metodo', 'modelo', 'n', 'parametros'] + \
            [c for x in nombres_x for c in (x, 's_' + x)] + \
            ['phi0', 'lambda0', 'media_l', 'std_l', 'sigma0', 'rmse_residuos', 'R2', 'R2_ajustado', 'mejor']
        escribir_csv(os.path.join(args.salida, f'modelo_parametrico_{m.nombre}.csv'), campos_res, filas_res)

        filas_int = []
        idx_com = {p['id']: k for k, p in enumerate(comunes)}
        for p in puntos:
            r = {'id': p['id'], 'lat': fmt(p['lat'], 9), 'lon': fmt(p['lon'], 9),
                 'h_ZT': fmt(p['h_zt']), 'H_ZT': fmt(p['H_zt']),
                 'usado_en_modelo': 'SI' if p['id'] in idx_com else 'NO'}
            for met in METODOS:
                c = CLAVE_METODO[met]
                z = res[met].get(p['id'])
                r['N_' + c] = fmt(z)
                r['l_' + c] = fmt(p['h_zt'] - p['H_zt'] - z) if z is not None else ''
                r['v_' + c] = fmt(ajustes[met]['v'][idx_com[p['id']]]) if p['id'] in idx_com else ''
            r['sigma_kriging'] = fmt(sig_k[p['id']]['sigma']) if p['id'] in sig_k else ''
            r['variograma_kriging'] = sig_k[p['id']]['modelo'] if p['id'] in sig_k else ''
            filas_int.append(r)
        campos_int = ['id', 'lat', 'lon', 'h_ZT', 'H_ZT', 'usado_en_modelo'] + \
            [f'{pre}_{CLAVE_METODO[met]}' for pre in ('N', 'l', 'v') for met in METODOS] + \
            ['sigma_kriging', 'variograma_kriging']
        escribir_csv(os.path.join(args.salida, f'interpolacion_{m.nombre}.csv'), campos_int, filas_int)

    campos_g = ['malla', 'metodo', 'modelo', 'n', 'parametros'] + \
        [c for x in (['dH', 'x1', 'x2'] if args.modelo == '18' else ['x0', 'x1', 'x2', 'x3'])
         for c in (x, 's_' + x)] + ['sigma0', 'rmse_residuos', 'R2', 'R2_ajustado', 'mejor']
    escribir_csv(os.path.join(args.salida, 'resumen_modelos.csv'), campos_g, resumen_global)
    log('\n' + '=' * 78)
    log(f'Archivos en {args.salida}: coordenadas_zero_tide.csv, puntos_rechazados.csv, '
        f'interpolacion_<malla>.csv, modelo_parametrico_<malla>.csv, resumen_modelos.csv, reporte.txt')
    log(f'Tiempo total: {time.time() - t0:.1f} s')
    with open(os.path.join(args.salida, 'reporte.txt'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(BITACORA) + '\n')


if __name__ == '__main__':
    main()
