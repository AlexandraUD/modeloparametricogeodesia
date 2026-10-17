#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interfaz gráfica: Puntos Red Nacional de Nivelación frente a modelos de geoide.

Solo usa Python estándar (tkinter). Los cálculos los hace
parametrico1.py (debe estar en la misma carpeta), así que la
interfaz y la línea de comandos dan exactamente los mismos resultados.

Ejecutar:   python parametrico2.py

Referencias
  Guimarães, G. N., Matos, A. C. O. C., Blitzkow, D. & Almeida, F. G. V. (2025).
  High-resolution regional gravity field modeling to connect the local vertical
  systems in South America to the IHRF. Earth, Planets and Space 77:98.
  https://doi.org/10.1186/s40623-025-02226-5
  EGM2008: ICGEM (GFZ) - https://icgem.gfz.de/calc_res/d533101e8204f4397361b8d33267c6ed173ceea602f93b06f39678a51579bc4c
  QGeoidCOL2023 y coords_QGEOID2023: IGAC - https://redgeodesica.igac.gov.co/
"""

import math
import os
import queue
import subprocess
import sys
import threading
import types
import webbrowser
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import parametrico1 as mp

AQUI = os.path.dirname(os.path.abspath(__file__))
URL_PAPER = 'https://doi.org/10.1186/s40623-025-02226-5'
URL_ICGEM = ('https://icgem.gfz.de/calc_res/'
             'd533101e8204f4397361b8d33267c6ed173ceea602f93b06f39678a51579bc4c')
URL_IGAC = 'https://redgeodesica.igac.gov.co/'
CITA = ('Guimarães, G. N., Matos, A. C. O. C., Blitzkow, D. & Almeida, F. G. V. (2025). '
        'High-resolution regional gravity field modeling to connect the local vertical systems '
        'in South America to the IHRF. Earth, Planets and Space, 77, 98.')

# Paleta
AZUL = '#1f4e79'
AZUL2 = '#2e75b6'
FONDO = '#f3f6fa'
BLANCO = '#ffffff'
GRIS = '#1b2631'          # textos de ayuda: casi negro (antes gris)
NEGRO = '#1b2631'
FONDO_FORMULA = '#eef4fb'
FUENTE_MATE = 'Cambria'   # fuente con buena tipografía matemática en Windows
NOMBRE_PUNTOS = 'Puntos Red Nacional de Nivelación'
NOMBRE_MODELO = {'18': 'Modelo de desplazamiento e inclinación',
                 '23': 'Superficie de corrección de 4 parámetros'}
VERDE = '#1e8449'
VERDE_CLARO = '#d5f5e3'
ROJO = '#c0392b'
COLORES_MODELO = {'sph': '#2e75b6', 'exp': '#e67e22', 'gau': '#8e44ad', 'cir': '#16a085'}

ELEGIR = '— elija —'
NO_VIENE = '(no viene)'
SIS_TEXTO = {'zero_tide': 'Zero tide', 'mean_tide': 'Mean tide', 'tide_free': 'Tide free'}
TEXTO_SIS = {v: k for k, v in SIS_TEXTO.items()}


# Las explicaciones usan una marcación mínima: $...$ = fórmula, _{..} = subíndice,
# ^{..} = superíndice (ver pintar()).
def explicar_h(s):
    return {None: 'Elija el sistema de mareas de h.',
            'zero_tide': 'Ya está en zero tide: no se modifica.',
            'tide_free': 'Se corrige así:   $h_{ZT} = h_{TF} + 0.0602 − 0.1790 sin²φ − 0.0018 sin⁴φ$',
            'mean_tide': 'En mean tide la corteza es la misma que en zero tide: no se modifica.'}[s]


def explicar_H(s):
    return {None: 'Elija el sistema de mareas de H.',
            'zero_tide': 'Ya está en zero tide: no se modifica.',
            'mean_tide': 'Se corrige así:   $H_{ZT} = H_{MT} − 0.198 ( ³⁄₂ sin²φ − ½ )$',
            'tide_free': 'Se corrige así:   $H_{ZT} = H_{TF} + Δh − ΔN$   (Δh y ΔN: correcciones de h y de N; '
                         'caso poco habitual)'}[s]


def explicar_N(s):
    return {None: 'Elija el sistema de mareas de la malla.',
            'zero_tide': 'Ya está en zero tide: no se modifica.',
            'tide_free': 'Se corrige así:   $N_{ZT} = N_{TF} − 0.0594 ( ³⁄₂ sin²φ − ½ )$   y se guarda una copia '
                         'de la malla corregida (<malla>_zero_tide.txt).',
            'mean_tide': 'Se corrige así:   $N_{ZT} = N_{MT} + 0.198 ( ³⁄₂ sin²φ − ½ )$   y se guarda una copia '
                         'de la malla corregida (<malla>_zero_tide.txt).'}[s]


def preparar_tags(txt, tam=11, tam_mate=13):
    """Estilos de texto para pintar(): fórmula, subíndice y superíndice."""
    txt.tag_configure('m', font=(FUENTE_MATE, tam_mate))
    txt.tag_configure('sub', font=(FUENTE_MATE, max(tam_mate - 4, 8)), offset=-4)
    txt.tag_configure('sup', font=(FUENTE_MATE, max(tam_mate - 4, 8)), offset=6)
    txt.tag_configure('subt', font=('Segoe UI', max(tam - 3, 8)), offset=-3)


def pintar(txt, s, *tags):
    """Inserta s en un Text interpretando $...$ (fórmula), _{..} y ^{..}."""
    i, mate, buf = 0, False, ''

    def volcar():
        if buf:
            txt.insert('end', buf, tags + (('m',) if mate else ()))
    while i < len(s):
        c = s[i]
        if c == '$':
            volcar()
            buf, mate = '', not mate
            i += 1
        elif c in '_^' and s[i + 1:i + 2] == '{':
            volcar()
            buf = ''
            fin = s.index('}', i + 2)
            estilo = ('sub' if c == '_' else 'sup') if mate else 'subt'
            txt.insert('end', s[i + 2:fin], tags + (('m',) if mate else ()) + (estilo,))
            i = fin + 1
        else:
            buf += c
            i += 1
    volcar()


class CajaFormula(tk.Text):
    """Recuadro de solo lectura para explicaciones con fórmulas bien escritas."""

    def __init__(self, padre, alto=2, ancho=70):
        super().__init__(padre, height=alto, width=ancho, wrap='word', bg=FONDO_FORMULA, fg=NEGRO,
                         relief='flat', padx=10, pady=6, font=('Segoe UI', 11), cursor='arrow',
                         highlightthickness=0)
        preparar_tags(self)
        self.configure(state='disabled')

    def mostrar(self, s):
        self.configure(state='normal')
        self.delete('1.0', 'end')
        pintar(self, s)
        self.configure(state='disabled')


def fuente_de(nombre):
    n = nombre.lower()
    if 'egm' in n:
        return 'ICGEM (GFZ Potsdam)'
    if 'qgeoid' in n or 'geocol' in n:
        return 'IGAC – Red Geodésica'
    return '—'


def _turbo(t):
    """Paleta 'turbo' (aproximación polinómica de Google, 2019); t en [0, 1]."""
    r = 0.13572138 + t * (4.61539260 + t * (-42.66032258 + t * (132.13108234 + t * (-152.94239396 + t * 59.28637943))))
    g = 0.09140261 + t * (2.19418839 + t * (4.84296658 + t * (-14.18503333 + t * (4.27729857 + t * 2.82956604))))
    b = 0.10667330 + t * (12.64194608 + t * (-60.58204836 + t * (110.36276771 + t * (-89.90310912 + t * 27.34824973))))
    c = lambda v: max(0, min(255, int(round(v * 255))))
    return f'#{c(r):02x}{c(g):02x}{c(b):02x}'


PALETA = [_turbo(k / 255) for k in range(256)]


def curvas_de_nivel(m, niveles):
    """Marching squares sobre la malla: {nivel: [(lon1, lat1, lon2, lat2), ...]}."""
    import bisect as _b
    segs = {L: [] for L in niveles}
    lat, lon, Z = m.lat, m.lon, m.Z
    for i in range(len(lat) - 1):
        f0, f1 = Z[i], Z[i + 1]
        y0, y1 = lat[i], lat[i + 1]
        for j in range(len(lon) - 1):
            a, b, c, d = f0[j], f0[j + 1], f1[j + 1], f1[j]       # SO, SE, NE, NO
            if a is None or b is None or c is None or d is None:
                continue
            lo = min(a, b, c, d)
            hi = max(a, b, c, d)
            k0 = _b.bisect_right(niveles, lo)
            k1 = _b.bisect_right(niveles, hi)
            if k0 == k1:
                continue
            x0, x1 = lon[j], lon[j + 1]
            for L in niveles[k0:k1]:
                pts = []
                for (za, xa, ya), (zb, xb, yb) in (((a, x0, y0), (b, x1, y0)), ((b, x1, y0), (c, x1, y1)),
                                                   ((c, x1, y1), (d, x0, y1)), ((d, x0, y1), (a, x0, y0))):
                    if (za >= L) != (zb >= L):
                        t = (L - za) / (zb - za)
                        pts.append((xa + t * (xb - xa), ya + t * (yb - ya)))
                if len(pts) >= 2:
                    segs[L].append((pts[0][0], pts[0][1], pts[1][0], pts[1][1]))
                if len(pts) == 4:
                    segs[L].append((pts[2][0], pts[2][1], pts[3][0], pts[3][1]))
    return segs


def numeros_bonitos(vmin, vmax, n=6):
    """Marcas 'redondas' para un eje."""
    if vmax <= vmin:
        vmax = vmin + 1.0
    paso = (vmax - vmin) / max(n - 1, 1)
    e = math.floor(math.log10(paso))
    f = paso / 10 ** e
    paso = (1 if f < 1.5 else 2 if f < 3 else 5 if f < 7 else 10) * 10 ** e
    ini = math.floor(vmin / paso) * paso
    marcas, x = [], ini
    while x <= vmax + paso * 1e-9:
        marcas.append(round(x, 12))
        x += paso
    return marcas


def fmt_num(x):
    if x == 0:
        return '0'
    if abs(x) >= 1000 or abs(x) < 0.001:
        return f'{x:.2e}'
    return f'{x:.4g}'


class App:
    def __init__(self, root):
        self.root = root
        root.title('Interpolación Puntos Red Nacional de Nivelación')
        root.geometry('1280x840')
        root.minsize(1050, 700)
        root.configure(bg=FONDO)
        self.q = queue.Queue()
        self.hilo = None
        self.detener = False
        self.encabezado = []
        self.mallas = []          # dicts: ruta, nombre, icgem, cabecera, sistema
        self.salida = None        # resultado de mp.ejecutar
        self.ruta_leida = None
        self._estilos()
        self._encabezado()
        self.nb = ttk.Notebook(root)
        self.nb.pack(fill='both', expand=True, padx=10, pady=(6, 0))
        self._tab_guia()
        self._tab_puntos()
        self._tab_mallas()
        self._tab_mapa()
        self._tab_calcular()
        self._tab_resultados()
        self._tab_kriging()
        self.estado = tk.StringVar(value='Listo. Lea la guía (pestaña 1) y siga los pasos en orden.')
        tk.Label(root, textvariable=self.estado, anchor='w', bg=AZUL, fg=BLANCO,
                 font=('Segoe UI', 10), padx=10, pady=4).pack(fill='x', side='bottom')
        self._precargar()
        self.root.after(100, self._leer_cola)

    # ------------------------------------------------------------------ estilo
    def _estilos(self):
        st = ttk.Style()
        try:
            st.theme_use('clam')
        except tk.TclError:
            pass
        st.configure('.', font=('Segoe UI', 10), background=FONDO)
        st.configure('TFrame', background=FONDO)
        st.configure('Card.TFrame', background=BLANCO)
        st.configure('TLabel', background=FONDO)
        st.configure('Card.TLabel', background=BLANCO)
        st.configure('Paso.TLabel', background=FONDO, foreground=AZUL, font=('Segoe UI', 13, 'bold'))
        st.configure('Ayuda.TLabel', background=FONDO, foreground=NEGRO, font=('Segoe UI', 11))
        st.configure('TLabel', background=FONDO, foreground=NEGRO, font=('Segoe UI', 11))
        st.configure('TCheckbutton', background=FONDO, font=('Segoe UI', 10))
        st.configure('TRadiobutton', background=FONDO, font=('Segoe UI', 10))
        st.configure('Formula.TLabel', background=FONDO_FORMULA, foreground=NEGRO,
                     font=('Consolas', 11), padding=6)
        st.configure('Mejor.TLabel', background=VERDE_CLARO, foreground=VERDE,
                     font=('Segoe UI', 11, 'bold'), padding=8)
        st.configure('TLabelframe', background=FONDO)
        st.configure('TLabelframe.Label', background=FONDO, foreground=AZUL,
                     font=('Segoe UI', 10, 'bold'))
        st.configure('TNotebook', background=FONDO)
        st.configure('TNotebook.Tab', padding=(14, 6), font=('Segoe UI', 10, 'bold'))
        st.map('TNotebook.Tab', background=[('selected', BLANCO)], foreground=[('selected', AZUL)])
        st.configure('Accion.TButton', font=('Segoe UI', 11, 'bold'), foreground=BLANCO,
                     background=AZUL2, padding=(16, 8))
        st.map('Accion.TButton', background=[('active', AZUL), ('disabled', '#a9b7c6')])
        st.configure('Treeview', rowheight=26, font=('Segoe UI', 10))
        st.configure('Treeview.Heading', font=('Segoe UI', 10, 'bold'), foreground=AZUL)
        st.configure('Mini.TButton', font=('Segoe UI', 10), padding=(6, 1))
        st.configure('green.Horizontal.TProgressbar', background=VERDE)

    def _encabezado(self):
        f = tk.Frame(self.root, bg=AZUL)
        f.pack(fill='x')
        tk.Label(f, text='Análisis de la Discrepancia Entre las Alturas de la Red Nacional de Nivelación y los Modelos de Geoide EGM2008 y de Cuasigeoide QGeoidCOL2023', bg=AZUL,
                 fg=BLANCO, font=('Segoe UI', 17, 'bold')).pack(anchor='w', padx=16, pady=(12, 2))
        tk.Label(f, text='Siga los pasos: cargue los puntos a interpolar y las mallas de los modelos de geoide. '
                         'El programa estandariza las alturas a sistema de mareas zero tide, interpola las ondulaciones por seis métodos '
                         'y le indica cuál se ajusta mejor.',
                 bg=AZUL, fg='#e3edf7', font=('Segoe UI', 11)).pack(anchor='w', padx=16, pady=(0, 12))

    def _paso(self, padre, titulo, ayuda):
        ttk.Label(padre, text=titulo, style='Paso.TLabel').pack(anchor='w', padx=14, pady=(12, 0))
        ttk.Label(padre, text=ayuda, style='Ayuda.TLabel', wraplength=1180,
                  justify='left').pack(anchor='w', padx=14, pady=(2, 8))

    # ------------------------------------------------------------- 1. Guía
    def _tab_guia(self):
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  1 · Guía y teoría  ')
        cont = ttk.Frame(fr)
        cont.pack(fill='both', expand=True, padx=12, pady=10)
        txt = tk.Text(cont, wrap='word', bg=BLANCO, fg=NEGRO, relief='flat', padx=36, pady=22,
                      font=('Segoe UI', 11), spacing1=2, spacing2=2, spacing3=4, cursor='arrow')
        sb = ttk.Scrollbar(cont, command=txt.yview)
        txt.configure(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y')
        txt.pack(side='left', fill='both', expand=True)
        txt.tag_configure('h1', font=('Segoe UI', 18, 'bold'), foreground=AZUL, spacing3=8)
        txt.tag_configure('h2', font=('Segoe UI', 14, 'bold'), foreground=AZUL, spacing1=22, spacing3=6)
        txt.tag_configure('intro', font=('Segoe UI', 12), spacing3=6)
        txt.tag_configure('b', font=('Segoe UI', 11, 'bold'))
        txt.tag_configure('li', lmargin1=26, lmargin2=44)
        txt.tag_configure('card', background='#f4f7fb', lmargin1=24, lmargin2=24, rmargin=24)
        txt.tag_configure('ct', font=('Segoe UI', 12, 'bold'), foreground=AZUL2, spacing1=10)
        txt.tag_configure('fb', justify='center', background=FONDO_FORMULA, spacing1=10, spacing3=10,
                          lmargin1=60, lmargin2=60, rmargin=60)
        txt.tag_configure('def', lmargin1=60, lmargin2=110)
        txt.tag_configure('ref', lmargin1=24, lmargin2=44, spacing3=8)
        preparar_tags(txt, 11, 15)

        def link(texto, url, *tags):
            tag = f'link{len(txt.tag_names())}'
            txt.tag_configure(tag, foreground=AZUL2, underline=True)
            txt.tag_bind(tag, '<Button-1>', lambda e: webbrowser.open(url))
            txt.tag_bind(tag, '<Enter>', lambda e: txt.configure(cursor='hand2'))
            txt.tag_bind(tag, '<Leave>', lambda e: txt.configure(cursor='arrow'))
            txt.insert('end', texto, tags + (tag,))

        P = lambda s, *t: pintar(txt, s + '\n', *t)
        H2 = lambda s: txt.insert('end', s + '\n', ('h2',))
        F = lambda s: pintar(txt, '$' + s + '$\n', 'fb')

        def LI(*items):
            for it in items:
                pintar(txt, '•   ' + it + '\n', 'li')

        def CARD(titulo, *parrafos):
            txt.insert('end', titulo + '\n', ('card', 'ct'))
            for p in parrafos:
                pintar(txt, p + '\n', 'card')

        def DEF(simbolo, texto):
            pintar(txt, f'${simbolo}$', 'def')
            pintar(txt, '   ' + texto + '\n', 'def')

        txt.insert('end', 'Guía rápida\n', ('h1',))
        P('Este programa compara las alturas de los ' + NOMBRE_PUNTOS + ' con modelos de geoide '
          'y responde dos preguntas: ¿cuánto está desplazado el nivel de referencia de las alturas '
          'respecto al geoide? y ¿qué método de interpolación del geoide ofrece el mejor ajuste?', 'intro')

        H2('Pasos')
        LI('2 · Puntos — cargue el archivo de puntos a interpolar e indique qué es cada columna.',
           '3 · Mallas — cargue los modelos del geoide de interés por ejemplo: QGeoidCOL2023, EGM2008, etc.',
           '4 · Mapa — compruebe la ubicación de los puntos y qué zonas cubre cada malla.',
           '5 · Calcular — revise los parámetros de las interpolaciones y presione «Ejecutar cálculo».',
           '6 · Resultados — observe el ajuste de cada método y cuál podría ser el mejor.',
           '7 · Kriging ordinario — explore el semivariograma usado por el kriging.')

        H2('Datos')
        CARD(NOMBRE_PUNTOS + '   (coords_QGEOID2023.txt · IGAC)',
             'Puntos medidos en campo. Cada uno tiene su posición (latitud φ y longitud λ) y dos alturas:',
             '      h  altura elipsoidal: altura sobre el elipsoide, la superficie matemática de referencia.',
             '      H  altura sobre el nivel del mar: la que se obtiene con la nivelación.')
        CARD('QGeoidCOL2023   (IGAC)',
             'Modelo del geoide (cuasigeoide) de Colombia calculado por el IGAC. Es una malla de cuadrícula '
             'regular de nodos cada 5′ (unos 9 km) que cubre el país y sus mares. En cada nodo da la ondulación '
             'N, es decir, la separación entre el elipsoide y el geoide.')
        CARD('EGM2008   (ICGEM)',
             'Modelo global del geoide (Earth Gravitational Model 2008), calculado con datos de gravedad '
             'terrestres y satelitales para todo el planeta. Esta malla se generó en el servicio ICGEM para la '
             'región de Colombia, con nodos cada 2′ (unos 3.7 km). Su encabezado indica que viene en tide free.')
        P('')
        P('Puntos y mallas no son lo mismo: los puntos son mediciones; las mallas son modelos. Como un punto '
          'casi nunca cae justo sobre un nodo, el valor de N en cada punto se interpola a partir de los nodos '
          'vecinos. El nombre de los archivos no importa: lo que cuenta es cargar los puntos en «2 · Puntos» y '
          'las mallas en «3 · Mallas».')

        H2('Discrepancia l')
        P('Si las alturas y el geoide fueran perfectos, se cumpliría h = H + N. En la práctica queda una '
          'diferencia, la discrepancia l:')
        F('l = h − H − N')
        P('Si l fuera igual en todos los puntos, el nivel de referencia de las alturas estaría solo desplazado. '
          'Si cambia de un lugar a otro, hay inclinaciones o errores en la red o en el modelo.')

        H2('Sistema de mareas')
        P('La atracción del Sol y de la Luna deforma la Tierra de manera permanente. Según cómo se trate ese '
          'efecto existen tres convenciones: tide free, mean tide y zero tide. Para comparar h, H y N deben '
          'estar en la misma; aquí se usa zero tide. Usted indica en qué sistema viene cada dato y el programa '
          'corrige solo lo que haga falta (lo que ya está en zero tide no se toca):')
        F('H_{ZT} = H_{MT} − 0.198 ( ³⁄₂ sin²φ − ½ )')
        F('h_{ZT} = h_{TF} + 0.0602 − 0.1790 sin²φ − 0.0018 sin⁴φ')
        F('N_{ZT} = N_{TF} − 0.0594 ( ³⁄₂ sin²φ − ½ )')
        P('Los subíndices indican el sistema: ZT zero tide, MT mean tide, TF tide free. Son correcciones de '
          'pocos centímetros; por ejemplo, EGM2008 sube entre 2.5 y 3 cm en Colombia al pasarlo a zero tide.')

        H2('Depuración de datos')
        P('Antes de calcular se revisa cada punto y se descartan los que tienen errores evidentes y se especifican en el archivo "puntos_rechazados.xlsx", dentro de la '
          'carpeta de resultados, con el motivo. También se ven en la pestaña «6 · Resultados». Se descarta un '
          'punto cuando:')
        LI('le falta un dato o un valor no es un número;',
           'su latitud o longitud son imposibles;',
           'está repetido (mismo identificador o mismas coordenadas);',
           'su altura elipsoidal vale cero (dato faltante);',
           'h − H es mayor de 150 m (ningún geoide se separa tanto);',
           'su discrepancia l se aleja demasiado de la del resto de puntos (error grueso: más de 5 veces la '
           'dispersión típica; el umbral se puede cambiar en «5 · Calcular»);',
           'queda fuera de la zona que cubre alguna malla.')

        H2('Interpolación con 6 métodos')
        P('Para obtener N en cada punto se prueban seis formas de interpolar la malla:')
        LI('Vecino más cercano — toma el valor del nodo más próximo.',
           'Lineal por triángulos — divide cada celda en dos triángulos y usa un plano.',
           'Bilineal — combina los 4 nodos de la celda según la cercanía.',
           'Bicúbica — usa 16 nodos y una curva suave.',
           'IDW (inverso de la distancia) — promedia los nodos cercanos dando más peso a los más próximos; '
           'el exponente p controla cuánto más. En «5 · Calcular» el programa puede sugerir el mejor p.',
           'Kriging ordinario — calcula los pesos a partir de cómo cambia N con la distancia (ver abajo).')

        H2('Kriging ordinario y semivariograma')
        P('El kriging primero estudia cómo varía N con la distancia. Para eso arma el semivariograma: compara '
          'pares de nodos y mide cuánto difieren según lo separados que estén; normalmente, cuanto más lejos, '
          'más diferentes. El covariograma muestra lo mismo al revés: cuánto se parecen.')
        P('A los puntos experimentales del semivariograma se les ajustan cuatro curvas (esférica, exponencial, '
          'gaussiana y circular) y se elige la de menor SSRw.')
        CARD('¿Qué es SSRw?',
             'Es la suma de cuadrados ponderada de los residuos: mide qué tan lejos queda la curva de los puntos '
             'experimentales. Cada diferencia se eleva al cuadrado y se multiplica por el número de pares de '
             'nodos de ese punto, porque un punto calculado con más pares es más confiable. Menor SSRw = la curva '
             'sigue mejor a los datos.')
        P('')
        P('Como el geoide cambia mucho de una región a otra, el semivariograma se ajusta en una ventana local '
          'alrededor de cada punto. Los puntos que caen en la misma celda de la malla comparten ventana, así '
          'que se ajusta un semivariograma por celda. En «7 · Kriging ordinario» puede ver un resumen de '
          'todos o el de un punto en particular.')

        H2('Modelo de ajuste')
        P('Con la discrepancia l de todos los puntos se ajusta, por mínimos cuadrados, una superficie sencilla '
          'que recoge su parte sistemática. Hay dos opciones:')
        CARD('Modelo de desplazamiento e inclinación   (recomendado)',
             'Supone que la diferencia es un desplazamiento vertical común más una inclinación suave en sentido '
             'norte–sur y otra en sentido este–oeste:')
        F('l = δH + x_{1} (φ − φ_{0}) + x_{2} (λ − λ_{0}) cos φ')
        DEF('δH', 'desplazamiento vertical medio (m): cuánto está corrido, en promedio, el nivel de referencia '
                  'de las alturas respecto al geoide.')
        DEF('x_{1}', 'inclinación norte–sur (m por grado de latitud).')
        DEF('x_{2}', 'inclinación este–oeste (m por grado de longitud).')
        DEF('φ_{0}, λ_{0}', 'latitud y longitud medias de los puntos (el centro de la red).')
        P('')
        CARD('Superficie de corrección de 4 parámetros',
             'Una superficie algo más flexible, habitual para ajustar un geoide a puntos nivelados:')
        F('l = x_{0} + x_{1} cos φ cos λ + x_{2} cos φ sin λ + x_{3} sin φ')
        DEF('x_{0} … x_{3}', 'coeficientes en metros que, juntos, describen la superficie.')
        P('')
        P('En «6 · Resultados» verá cada parámetro con su desviación estándar (±), σ₀ (cuánto se separan en '
          'promedio los puntos de la superficie ajustada, en metros), R² y R² ajustado.')

        H2('R² ajustado')
        P('R² indica qué parte de la variación de la discrepancia l explica el modelo: 0 es nada y 1 es todo. '
          'El R² ajustado lo corrige según el número de parámetros, para no favorecer a los modelos con más '
          'coeficientes:')
        F('R²_{aj} = 1 − (1 − R²) · (n − 1) ⁄ (n − u)')
        P('n es el número de puntos y u el número de parámetros (3 o 4). El modelo se ajusta seis veces, una '
          'con el N de cada método y siempre con los mismos puntos: el método con mayor R² ajustado es el mejor.')

        H2('Archivos que se guardan')
        LI('coordenadas_zero_tide.xlsx — los puntos aceptados con sus alturas en zero tide.',
           'puntos_rechazados.xlsx — los puntos descartados y el motivo.',
           'resumen_modelos.xlsx — parámetros, σ₀, R² y R² ajustado de cada método y malla.',
           'interpolacion_<malla>.csv — N, l y residuos de los 6 métodos en cada punto.',
           'reporte.txt — la bitácora completa del cálculo.')

        H2('Referencias')
        P('Las correcciones de marea y los modelos de ajuste siguen a:')
        txt.insert('end', CITA + '  ', ('ref',))
        link(URL_PAPER, URL_PAPER, 'ref')
        txt.insert('end', '\n')
        txt.insert('end', 'EGM2008: Pavlis, N. K., Holmes, S. A., Kenyon, S. C. & Factor, J. K. (2012). The development '
                          'and evaluation of the Earth Gravitational Model 2008 (EGM2008). Journal of Geophysical '
                          'Research, 117, B04406. Malla generada en ICGEM (GFZ Potsdam): ', ('ref',))
        link(URL_ICGEM, URL_ICGEM, 'ref')
        txt.insert('end', '\n')
        txt.insert('end', 'QGeoidCOL2023 y ' + NOMBRE_PUNTOS + ': Instituto Geográfico Agustín Codazzi '
                          '(IGAC), Red Geodésica Nacional: ', ('ref',))
        link(URL_IGAC, URL_IGAC, 'ref')
        txt.insert('end', '\n')
        txt.configure(state='disabled')

    # ------------------------------------------------------------ 2. Puntos
    def _tab_puntos(self):
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  2 · Puntos  ')
        self._paso(fr, 'Paso 1 — ' + NOMBRE_PUNTOS,
                   'Son los puntos medidos (por ejemplo, coords_QGEOID2023.txt del IGAC). En ellos se interpola N '
                   'y se calcula la discrepancia l. Formatos: texto con encabezado separado por «;», «,», '
                   'tabulación o espacios.')
        f1 = ttk.Frame(fr)
        f1.pack(fill='x', padx=14)
        self.v_puntos = tk.StringVar()
        ttk.Entry(f1, textvariable=self.v_puntos).pack(side='left', fill='x', expand=True)
        ttk.Button(f1, text='Examinar…', command=self._elegir_puntos).pack(side='left', padx=6)
        ttk.Button(f1, text='Leer', command=lambda: self._cargar_puntos(self.v_puntos.get())).pack(side='left')

        self.lbl_info_puntos = ttk.Label(fr, text='', style='Ayuda.TLabel')
        self.lbl_info_puntos.pack(anchor='w', padx=14, pady=(4, 0))
        vista = ttk.LabelFrame(fr, text='Vista previa (primeras filas)')
        vista.pack(fill='both', expand=True, padx=14, pady=6)
        self.tv_prev = ttk.Treeview(vista, show='headings', height=7)
        sx = ttk.Scrollbar(vista, orient='horizontal', command=self.tv_prev.xview)
        self.tv_prev.configure(xscrollcommand=sx.set)
        self.tv_prev.pack(fill='both', expand=True, padx=4, pady=(4, 0))
        sx.pack(fill='x', padx=4)

        abajo = ttk.Frame(fr)
        abajo.pack(fill='x', padx=14, pady=(0, 10))
        col = ttk.LabelFrame(abajo, text='¿Qué es cada columna?')
        col.pack(side='left', fill='both', expand=True, padx=(0, 8))
        self.cb_col = {}
        textos = [('id', 'Identificador (opcional)'), ('lat', 'Latitud φ [grados]'),
                  ('lon', 'Longitud λ [grados]'), ('h', 'Altura elipsoidal h [m]'),
                  ('H', 'Altura sobre el nivel del mar H [m] (nivelación)')]
        for r, (k, t) in enumerate(textos):
            ttk.Label(col, text=t).grid(row=r, column=0, sticky='w', padx=8, pady=3)
            cb = ttk.Combobox(col, state='readonly', width=28)
            cb.grid(row=r, column=1, sticky='w', padx=8, pady=3)
            self.cb_col[k] = cb
        ttk.Label(col, text='Para el análisis se necesitan h y H.', style='Ayuda.TLabel', wraplength=420
                  ).grid(row=5, column=0, columnspan=2, sticky='w', padx=8, pady=(2, 6))

        mar = ttk.LabelFrame(abajo, text='¿En qué sistema de mareas vienen las alturas?')
        mar.pack(side='left', fill='both', expand=True)
        self.cb_sis_h = ttk.Combobox(mar, state='readonly', width=16, values=list(SIS_TEXTO.values()))
        self.cb_sis_H = ttk.Combobox(mar, state='readonly', width=16, values=list(SIS_TEXTO.values()))
        self.lbl_h = CajaFormula(mar, alto=2, ancho=64)
        self.lbl_H = CajaFormula(mar, alto=2, ancho=64)
        self.lbl_h.mostrar(explicar_h(None))
        self.lbl_H.mostrar(explicar_H(None))
        ttk.Label(mar, text='h (elipsoidal):').grid(row=0, column=0, sticky='w', padx=8, pady=(6, 2))
        self.cb_sis_h.grid(row=0, column=1, sticky='w', pady=(6, 2))
        self.lbl_h.grid(row=1, column=0, columnspan=2, sticky='we', padx=8)
        ttk.Label(mar, text='H (nivel del mar):').grid(row=2, column=0, sticky='w', padx=8, pady=(10, 2))
        self.cb_sis_H.grid(row=2, column=1, sticky='w', pady=(10, 2))
        self.lbl_H.grid(row=3, column=0, columnspan=2, sticky='we', padx=8)
        ttk.Label(mar, text='Si no lo sabe, revise los metadatos de la fuente (IGAC). Lo usual es '
                            'h en tide free y H en mean tide.', style='Ayuda.TLabel', wraplength=560
                  ).grid(row=4, column=0, columnspan=2, sticky='w', padx=8, pady=6)
        self.cb_sis_h.bind('<<ComboboxSelected>>',
                           lambda e: self.lbl_h.mostrar(explicar_h(self._sis(self.cb_sis_h))))
        self.cb_sis_H.bind('<<ComboboxSelected>>',
                           lambda e: self.lbl_H.mostrar(explicar_H(self._sis(self.cb_sis_H))))

    @staticmethod
    def _sis(cb):
        return TEXTO_SIS.get(cb.get())

    def _elegir_puntos(self):
        r = filedialog.askopenfilename(title='Archivo de ' + NOMBRE_PUNTOS, initialdir=AQUI,
                                       filetypes=[('Texto / CSV', '*.txt *.csv *.dat'), ('Todos', '*.*')])
        if r:
            self.v_puntos.set(r)
            self._cargar_puntos(r)

    def _cargar_puntos(self, ruta):
        if not ruta:
            return
        try:
            enc, filas, delim = mp.leer_tabla(ruta)
        except SystemExit as e:
            messagebox.showerror('Archivo de puntos', str(e))
            return
        except Exception as e:
            messagebox.showerror('Archivo de puntos', f'No se pudo leer el archivo:\n{e}')
            return
        self.encabezado = enc
        self.ruta_leida = os.path.abspath(ruta)
        self.tv_prev.delete(*self.tv_prev.get_children())
        self.tv_prev.configure(columns=[f'c{k}' for k in range(len(enc))])
        for k, c in enumerate(enc):
            self.tv_prev.heading(f'c{k}', text=f'{k + 1}: {c}')
            self.tv_prev.column(f'c{k}', width=150, anchor='w', stretch=True)
        for fila in filas[:12]:
            self.tv_prev.insert('', 'end', values=fila)
        valores = [NO_VIENE] + enc
        for k, cb in self.cb_col.items():
            cb.configure(values=valores)
            sug = mp.sugerir_columna(enc, k)
            cb.set(sug if sug else NO_VIENE)
        sep = {';': '«;»', ',': '«,»', '\t': 'tabulación', None: 'espacios'}.get(delim, repr(delim))
        self.lbl_info_puntos.configure(
            text=f'{len(filas)} registros, {len(enc)} columnas, separador {sep}. Las columnas se '
                 f'sugirieron por su nombre: revíselas abajo.')
        self.estado.set(f'Puntos leídos: {os.path.basename(ruta)} ({len(filas)} registros).')

    # ------------------------------------------------------------ 3. Mallas
    def _tab_mallas(self):
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  3 · Mallas  ')
        self._paso(fr, 'Paso 2 — Mallas del geoide (N)',
                   'La malla es un modelo del geoide, es una cuadrícula regular de nodos con el valor de N, más no son '
                   'puntos medidos. De la malla se interpola N en los ' + NOMBRE_PUNTOS + ' con 6 métodos, y '
                   'cada malla se analiza por separado. Formatos: «lat lon N» sin encabezado (como QGeoidCOL2023 '
                   'del IGAC) o archivo de ICGEM con encabezado (como EGM2008).')
        caja = ttk.Frame(fr)
        caja.pack(fill='both', expand=True, padx=14)
        cols = ('archivo', 'formato', 'cabecera', 'sistema', 'fuente')
        self.tv_mallas = ttk.Treeview(caja, columns=cols, show='headings', height=6, selectmode='browse')
        for c, t, w in zip(cols, ('Archivo', 'Formato', 'Encabezado dice', 'Sistema de mareas elegido',
                                  'Fuente'), (330, 150, 150, 200, 200)):
            self.tv_mallas.heading(c, text=t)
            self.tv_mallas.column(c, width=w, anchor='w')
        self.tv_mallas.pack(fill='x')
        self.tv_mallas.tag_configure('falta', foreground=ROJO)
        self.tv_mallas.bind('<<TreeviewSelect>>', self._malla_sel)
        bot = ttk.Frame(caja)
        bot.pack(fill='x', pady=6)
        ttk.Button(bot, text='Agregar malla…', command=self._agregar_malla).pack(side='left')
        ttk.Button(bot, text='Quitar seleccionada', command=self._quitar_malla).pack(side='left', padx=6)

        sel = ttk.LabelFrame(caja, text='Sistema de mareas de la malla seleccionada')
        sel.pack(fill='x', pady=6)
        self.lbl_malla_sel = ttk.Label(sel, text='Seleccione una malla en la tabla.')
        self.lbl_malla_sel.grid(row=0, column=0, sticky='w', padx=8, pady=6)
        self.cb_sis_malla = ttk.Combobox(sel, state='readonly', width=16, values=list(SIS_TEXTO.values()))
        self.cb_sis_malla.grid(row=0, column=1, sticky='w', padx=8)
        self.cb_sis_malla.bind('<<ComboboxSelected>>', self._cambiar_sis_malla)
        self.lbl_N = CajaFormula(sel, alto=2, ancho=120)
        self.lbl_N.mostrar(explicar_N(None))
        self.lbl_N.grid(row=1, column=0, columnspan=2, sticky='we', padx=8, pady=(0, 8))
        ttk.Label(caja, text='El sistema de mareas lo elige el usuario.',
                  style='Ayuda.TLabel').pack(anchor='w', pady=4)

    def _agregar_malla(self, ruta=None, sistema=None):
        if ruta is None:
            ruta = filedialog.askopenfilename(title='Malla de N', initialdir=AQUI,
                                              filetypes=[('Mallas', '*.txt *.gdf *.dat *.grd'),
                                                         ('Todos', '*.*')])
            if not ruta:
                return
        if any(os.path.abspath(m['ruta']) == os.path.abspath(ruta) for m in self.mallas):
            messagebox.showinfo('Mallas', 'Esa malla ya está en la lista.')
            return
        try:
            icgem, cab = mp.leer_cabecera_malla(ruta)
        except OSError as e:
            messagebox.showerror('Mallas', f'No se pudo abrir la malla:\n{e}')
            return
        nombre = os.path.splitext(os.path.basename(ruta))[0]
        self.mallas.append({'ruta': ruta, 'nombre': nombre, 'icgem': icgem, 'cabecera': cab,
                            'sistema': sistema if sistema else cab})
        self._refrescar_mallas()

    def _quitar_malla(self):
        s = self.tv_mallas.selection()
        if s:
            del self.mallas[int(s[0])]
            self._refrescar_mallas()

    def _refrescar_mallas(self):
        self.tv_mallas.delete(*self.tv_mallas.get_children())
        for k, m in enumerate(self.mallas):
            sis = SIS_TEXTO.get(m['sistema'], ELEGIR)
            self.tv_mallas.insert('', 'end', iid=str(k), tags=() if m['sistema'] else ('falta',),
                                  values=(os.path.basename(m['ruta']),
                                          'ICGEM (lon, lat, N)' if m['icgem'] else 'lat lon N',
                                          SIS_TEXTO.get(m['cabecera'], '(no dice)'), sis,
                                          fuente_de(m['nombre'])))

    def _malla_sel(self, _=None):
        s = self.tv_mallas.selection()
        if not s:
            return
        m = self.mallas[int(s[0])]
        self.lbl_malla_sel.configure(text=f'{os.path.basename(m["ruta"])}:')
        self.cb_sis_malla.set(SIS_TEXTO.get(m['sistema'], ''))
        self.lbl_N.mostrar(explicar_N(m['sistema']))

    def _cambiar_sis_malla(self, _=None):
        s = self.tv_mallas.selection()
        if not s:
            messagebox.showinfo('Mallas', 'Primero seleccione una malla en la tabla.')
            return
        k = int(s[0])
        self.mallas[k]['sistema'] = self._sis(self.cb_sis_malla)
        self._refrescar_mallas()
        self.tv_mallas.selection_set(str(k))
        self.lbl_N.mostrar(explicar_N(self.mallas[k]['sistema']))

    # -------------------------------------------------------------- 4. Mapa
    def _tab_mapa(self):
        self.mapa_mallas = {}         # ruta -> Malla (valores tal como vienen en el archivo)
        self.mapa_cargando = set()
        self.mapa_curvas = {}         # (ruta, paso) -> segmentos
        self.mapa_img = None
        self.mapa_geo = None          # transformación pantalla <-> lon/lat
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  4 · Mapa  ')
        self._paso(fr, 'Mapa — dónde están los puntos y qué zona cubre cada malla',
                   'Rectángulos punteados: zona que cubre cada malla. Puntos negros: ' + NOMBRE_PUNTOS + ' (después del cálculo, los apartados en '
                   'la depuración aparecen como × rojas). Pase el ratón sobre el mapa para leer lat, lon y N.')
        ctrl = ttk.Frame(fr)
        ctrl.pack(fill='x', padx=14)
        ttk.Label(ctrl, text='Fondo (malla):').pack(side='left')
        self.cb_mapa_fondo = ttk.Combobox(ctrl, state='readonly', width=24)
        self.cb_mapa_fondo.pack(side='left', padx=(4, 12))
        self.cb_mapa_fondo.bind('<<ComboboxSelected>>', lambda e: self._dibujar_mapa())
        self.v_m = {k: tk.BooleanVar(value=v) for k, v in
                    (('puntos', True), ('cobertura', True), ('curvas', True), ('nodos', False))}
        for k, t in (('puntos', NOMBRE_PUNTOS), ('cobertura', 'Cobertura de las mallas'),
                     ('curvas', 'Curvas de nivel'), ('nodos', 'Nodos de la malla')):
            ttk.Checkbutton(ctrl, text=t, variable=self.v_m[k],
                            command=self._dibujar_mapa).pack(side='left', padx=4)
        ttk.Label(ctrl, text='  cada').pack(side='left')
        self.cb_mapa_paso = ttk.Combobox(ctrl, state='readonly', width=4, values=['2', '5', '10'])
        self.cb_mapa_paso.set('5')
        self.cb_mapa_paso.pack(side='left', padx=2)
        self.cb_mapa_paso.bind('<<ComboboxSelected>>', lambda e: self._dibujar_mapa())
        ttk.Label(ctrl, text='m').pack(side='left')
        ttk.Button(ctrl, text='Actualizar', command=self._dibujar_mapa).pack(side='right')

        cuerpo = ttk.Frame(fr)
        cuerpo.pack(fill='both', expand=True, padx=14, pady=8)
        self.cv_mapa = tk.Canvas(cuerpo, bg=BLANCO, highlightthickness=1, highlightbackground='#d5dbe3')
        self.cv_mapa.pack(side='left', fill='both', expand=True)
        self.cv_mapa.bind('<Configure>', lambda e: self._programar_mapa())
        self.cv_mapa.bind('<Motion>', self._mapa_cursor)
        self.nb.bind('<<NotebookTabChanged>>', lambda e: self._programar_mapa(), add='+')
        self._mapa_pendiente = None
        der = ttk.Frame(cuerpo)
        der.pack(side='left', fill='y', padx=(10, 0))
        self.lbl_cursor = ttk.Label(der, text='lat —   lon —   N —', style='Formula.TLabel', width=40)
        self.lbl_cursor.pack(fill='x')
        self.txt_mapa = tk.Text(der, width=46, font=('Consolas', 10), bg=BLANCO, relief='flat',
                                padx=10, pady=8, wrap='word')
        self.txt_mapa.pack(fill='both', expand=True, pady=(8, 0))

    def _programar_mapa(self):
        """Redibuja el mapa 200 ms después del último evento, solo si su pestaña está visible."""
        if self._mapa_pendiente is not None:
            self.root.after_cancel(self._mapa_pendiente)
        self._mapa_pendiente = self.root.after(200, self._mapa_si_visible)

    def _mapa_si_visible(self):
        self._mapa_pendiente = None
        if self.nb.index(self.nb.select()) == 3:
            self._dibujar_mapa()

    def _leer_coords_puntos(self):
        """Lat/lon del archivo de puntos con las columnas elegidas en el paso 1."""
        ruta = self.v_puntos.get().strip()
        if not ruta or not os.path.isfile(ruta) or not self.encabezado:
            return [], 'Sin archivo de puntos (paso 1).'
        ilat = self.cb_col['lat'].get()
        ilon = self.cb_col['lon'].get()
        iid = self.cb_col['id'].get()
        if ilat in ('', NO_VIENE) or ilon in ('', NO_VIENE):
            return [], 'Elija las columnas de latitud y longitud (paso 1).'
        clave = (ruta, ilat, ilon, iid)
        if getattr(self, '_cache_coords', (None,))[0] == clave:
            return self._cache_coords[1], None
        _, filas, _ = mp.leer_tabla(ruta)
        a, b = self.encabezado.index(ilat), self.encabezado.index(ilon)
        c = self.encabezado.index(iid) if iid not in ('', NO_VIENE) else None
        pts = []
        for k, f in enumerate(filas):
            la = mp._num(f[a]) if a < len(f) else None
            lo = mp._num(f[b]) if b < len(f) else None
            if la is None or lo is None:
                continue
            if lo > 180:
                lo -= 360
            if abs(la) <= 90 and abs(lo) <= 180:
                pts.append((f[c].strip() if c is not None and c < len(f) else str(k + 1), la, lo))
        self._cache_coords = (clave, pts)
        return pts, None

    def _malla_mapa(self, ruta):
        """Devuelve la malla ya leída o lanza su lectura en segundo plano."""
        if ruta in self.mapa_mallas:
            return self.mapa_mallas[ruta]
        if ruta not in self.mapa_cargando:
            self.mapa_cargando.add(ruta)
            self.estado.set(f'Leyendo {os.path.basename(ruta)} para el mapa…')

            def leer():
                try:
                    self.q.put(('mapa_ok', ruta, mp.cargar_malla(ruta)))
                except SystemExit as e:
                    self.q.put(('mapa_error', ruta, str(e)))
                except Exception as e:
                    self.q.put(('mapa_error', ruta, f'No se pudo leer la malla para el mapa:\n{e}'))
            threading.Thread(target=leer, daemon=True).start()
        return None

    def _dibujar_mapa(self):
        if not hasattr(self, 'cv_mapa'):
            return
        cv = self.cv_mapa
        W, H = cv.winfo_width(), cv.winfo_height()
        if W < 150 or H < 150:
            return
        # opciones de fondo
        nombres = [m['nombre'] for m in self.mallas]
        self.cb_mapa_fondo.configure(values=nombres + ['(ninguno)'])
        if self.cb_mapa_fondo.get() not in nombres + ['(ninguno)']:
            self.cb_mapa_fondo.set(nombres[0] if nombres else '(ninguno)')
        fondo = next((m for m in self.mallas if m['nombre'] == self.cb_mapa_fondo.get()), None)
        cargadas = {m['ruta']: self._malla_mapa(m['ruta']) for m in self.mallas}
        mf = cargadas.get(fondo['ruta']) if fondo else None
        pts, aviso_pts = self._leer_coords_puntos() if self.v_m['puntos'].get() else ([], None)

        cv.delete('all')
        # extensión: fondo + puntos (+ coberturas si se muestran)
        cajas = []
        if mf is not None:
            cajas.append((mf.lon[0], mf.lon[-1], mf.lat[0], mf.lat[-1]))
        if pts:
            cajas.append((min(p[2] for p in pts), max(p[2] for p in pts),
                          min(p[1] for p in pts), max(p[1] for p in pts)))
        if mf is None and (self.v_m['cobertura'].get() or not cajas):
            # sin fondo, el encuadre abarca todas las mallas; con fondo, el de esa malla
            cajas += [(m.lon[0], m.lon[-1], m.lat[0], m.lat[-1]) for m in cargadas.values() if m]
        if not cajas:
            msg = 'Cargue puntos (paso 1) o mallas (paso 2) para ver el mapa.'
            if self.mapa_cargando:
                msg = 'Leyendo las mallas… (unos segundos)'
            cv.create_text(W / 2, H / 2, text=msg, fill=GRIS, font=('Segoe UI', 11))
            return
        lon0, lon1 = min(c[0] for c in cajas), max(c[1] for c in cajas)
        lat0, lat1 = min(c[2] for c in cajas), max(c[3] for c in cajas)
        mx, my = (lon1 - lon0) * 0.02 or 0.05, (lat1 - lat0) * 0.02 or 0.05
        lon0, lon1, lat0, lat1 = lon0 - mx, lon1 + mx, lat0 - my, lat1 + my
        # márgenes y escala con aspecto geográfico correcto (1° lon = cos φ · 1° lat)
        L, R, T, B = 70, 110 if mf is not None else 20, 40, 52
        aw, ah = W - L - R, H - T - B
        k = math.cos(math.radians((lat0 + lat1) / 2))
        esc = min(aw / ((lon1 - lon0) * k), ah / (lat1 - lat0))
        pw, ph = (lon1 - lon0) * k * esc, (lat1 - lat0) * esc
        ox, oy = L + (aw - pw) / 2, T + (ah - ph) / 2
        X = lambda lo: ox + (lo - lon0) * k * esc
        Y = lambda la: oy + (lat1 - la) * esc
        self.mapa_geo = (ox, oy, pw, ph, lon0, lon1, lat0, lat1, mf)

        info = []
        # fondo: imagen de N (vecino más cercano por píxel)
        if mf is not None:
            vals = [z for fila in mf.Z for z in fila if z is not None]
            zmin, zmax = min(vals), max(vals)
            if zmax - zmin < 1e-12:
                zmax = zmin + 1.0
            iw, ih = max(int(pw), 1), max(int(ph), 1)
            img = tk.PhotoImage(width=iw, height=ih)
            dl, dlo = mf.lat[1] - mf.lat[0], mf.lon[1] - mf.lon[0]
            nl, nlo = len(mf.lat), len(mf.lon)
            cols_j = []
            for px in range(iw):
                lo = lon0 + (px + 0.5) / iw * (lon1 - lon0)
                j = round((lo - mf.lon[0]) / dlo)
                cols_j.append(j if 0 <= j < nlo else None)
            filas = []
            for py in range(ih):
                la = lat1 - (py + 0.5) / ih * (lat1 - lat0)
                i = round((la - mf.lat[0]) / dl)
                if not 0 <= i < nl:
                    filas.append('{' + ' '.join([BLANCO] * iw) + '}')
                    continue
                fz = mf.Z[i]
                fila = []
                for j in cols_j:
                    z = fz[j] if j is not None else None
                    fila.append(BLANCO if z is None else PALETA[int((z - zmin) / (zmax - zmin) * 255)])
                filas.append('{' + ' '.join(fila) + '}')
            img.put(' '.join(filas), to=(0, 0))
            self.mapa_img = img
            cv.create_image(ox, oy, image=img, anchor='nw')
            # curvas de nivel
            if self.v_m['curvas'].get():
                paso = float(self.cb_mapa_paso.get())
                clave = (fondo['ruta'], paso)
                if clave not in self.mapa_curvas:
                    niveles = [paso * n for n in range(math.ceil(zmin / paso), math.floor(zmax / paso) + 1)]
                    self.mapa_curvas[clave] = curvas_de_nivel(mf, niveles)
                for nivel, segs in self.mapa_curvas[clave].items():
                    ancho = 1.2 if abs(nivel % (2 * paso)) < 1e-9 else 0.6
                    for x1, y1, x2, y2 in segs:
                        cv.create_line(X(x1), Y(y1), X(x2), Y(y2), fill='#222222', width=ancho)
                    if segs:                         # etiqueta en el segmento central
                        x1, y1, x2, y2 = segs[len(segs) // 2]
                        cv.create_text(X((x1 + x2) / 2), Y((y1 + y2) / 2), text=f'{nivel:g}',
                                       font=('Segoe UI', 8), fill='#111111')
            # nodos (submuestreados si son muchos)
            if self.v_m['nodos'].get():
                s = max(1, math.ceil(math.sqrt(nl * nlo / 30000)))
                for i in range(0, nl, s):
                    yy = Y(mf.lat[i])
                    for j in range(0, nlo, s):
                        if mf.Z[i][j] is not None:
                            xx = X(mf.lon[j])
                            cv.create_rectangle(xx, yy, xx + 1, yy + 1, outline='#333333')
                info.append(f'Nodos dibujados: 1 de cada {s} en cada dirección' if s > 1 else
                            'Nodos dibujados: todos')
            # barra de color
            bx, by0, by1 = ox + pw + 22, oy, oy + ph
            for n in range(int(by1 - by0)):
                cv.create_line(bx, by1 - n, bx + 16, by1 - n,
                               fill=PALETA[int(n / max(by1 - by0 - 1, 1) * 255)])
            cv.create_rectangle(bx, by0, bx + 16, by1, outline='#555555')
            for v in numeros_bonitos(zmin, zmax, 7):
                if zmin <= v <= zmax:
                    yy = by1 - (v - zmin) / (zmax - zmin) * (by1 - by0)
                    cv.create_line(bx + 16, yy, bx + 20, yy, fill='#555555')
                    cv.create_text(bx + 23, yy, text=f'{v:g}', anchor='w', font=('Segoe UI', 9))
            cv.create_text(bx + 8, by0 - 12, text='N [m]', font=('Segoe UI', 10, 'bold'))
            paso_min = (mf.lat[1] - mf.lat[0]) * 60
            titulo = (f'{os.path.basename(fondo["ruta"])}: {nl} x {nlo} nodos, paso {paso_min:.4g}\'')
        else:
            cv.create_rectangle(ox, oy, ox + pw, oy + ph, fill='#f8fafc', outline='')
            titulo = 'Ubicación de los datos cargados'
        # cobertura de cada malla
        if self.v_m['cobertura'].get():
            colores = ['#c0392b', '#1a5276', '#7d3c98', '#117864']
            for n, m in enumerate(self.mallas):
                g = cargadas.get(m['ruta'])
                if g is None:
                    continue
                col = colores[n % len(colores)]
                x1, x2 = max(X(g.lon[0]), ox), min(X(g.lon[-1]), ox + pw)     # recorte al marco
                y1, y2 = max(Y(g.lat[-1]), oy), min(Y(g.lat[0]), oy + ph)
                if x2 <= x1 or y2 <= y1:
                    continue
                cv.create_rectangle(x1, y1, x2, y2, outline=col, width=2, dash=(8, 4))
                sale = (X(g.lon[0]) < ox - 1 or X(g.lon[-1]) > ox + pw + 1 or
                        Y(g.lat[-1]) < oy - 1 or Y(g.lat[0]) > oy + ph + 1)
                cv.create_text(x1 + 4, y1 + 4 + 14 * n, anchor='nw', fill=col,
                               font=('Segoe UI', 9, 'bold'),
                               text=os.path.basename(m['ruta']) + ('  (continúa fuera del mapa)' if sale else ''))
        # puntos
        rech = set()
        if self.salida:
            rech = {r['id'] for r in self.salida['rechazos']}
        n_rech = 0
        for pid, la, lo in pts:
            xx, yy = X(lo), Y(la)
            if pid in rech:
                n_rech += 1
                cv.create_line(xx - 3, yy - 3, xx + 4, yy + 4, fill=ROJO, width=2)
                cv.create_line(xx - 3, yy + 3, xx + 4, yy - 4, fill=ROJO, width=2)
            else:
                cv.create_oval(xx - 1.6, yy - 1.6, xx + 1.6, yy + 1.6, fill='#000000', outline='')
        # ejes
        cv.create_rectangle(ox, oy, ox + pw, oy + ph, outline='#555555')
        for v in numeros_bonitos(lon0, lon1, 8):
            if lon0 <= v <= lon1:
                cv.create_line(X(v), oy + ph, X(v), oy + ph + 5, fill='#555555')
                cv.create_text(X(v), oy + ph + 15, text=f'{v:g}', font=('Segoe UI', 9))
        for v in numeros_bonitos(lat0, lat1, 9):
            if lat0 <= v <= lat1:
                cv.create_line(ox - 5, Y(v), ox, Y(v), fill='#555555')
                cv.create_text(ox - 8, Y(v), text=f'{v:g}', anchor='e', font=('Segoe UI', 9))
        cv.create_text(ox + pw / 2, oy + ph + 34, text='Longitud (grados)', font=('Segoe UI', 10))
        cv.create_text(ox - 48, oy + ph / 2, text='Latitud (grados)', angle=90, font=('Segoe UI', 10))
        cv.create_text(ox + pw / 2, oy - 18, text=titulo, font=('Segoe UI', 11, 'bold'), fill='#1b2631')

        # panel de información
        t = self.txt_mapa
        t.configure(state='normal')
        t.delete('1.0', 'end')
        t.tag_configure('t', font=('Segoe UI', 10, 'bold'), foreground=AZUL)
        t.insert('end', NOMBRE_PUNTOS + '\n', 't')
        if pts:
            t.insert('end', f'{os.path.basename(self.v_puntos.get())}\n'
                            f'{len(pts)} puntos con coordenadas válidas\n'
                            f'lat [{min(p[1] for p in pts):.3f}, {max(p[1] for p in pts):.3f}]\n'
                            f'lon [{min(p[2] for p in pts):.3f}, {max(p[2] for p in pts):.3f}]\n')
            if self.salida:
                t.insert('end', f'Rechazados en la depuración: {n_rech} (× rojas)\n')
        else:
            t.insert('end', (aviso_pts or 'Puntos ocultos.') + '\n')
        t.insert('end', '\nMallas\n', 't')
        for m in self.mallas:
            g = cargadas.get(m['ruta'])
            t.insert('end', f'{os.path.basename(m["ruta"])}  ({fuente_de(m["nombre"])})\n')
            if g is None:
                t.insert('end', '  leyendo…\n')
                continue
            vals = [z for fila in g.Z for z in fila if z is not None]
            t.insert('end', f'  {len(g.lat)} x {len(g.lon)} nodos, paso {(g.lat[1] - g.lat[0]) * 60:.4g}\'\n'
                            f'  lat [{g.lat[0]:.3f}, {g.lat[-1]:.3f}]\n'
                            f'  lon [{g.lon[0]:.3f}, {g.lon[-1]:.3f}]\n'
                            f'  N [{min(vals):.3f}, {max(vals):.3f}] m\n')
            if pts:
                fuera = sum(1 for _, la, lo in pts
                            if not (g.lat[0] <= la <= g.lat[-1] and g.lon[0] <= lo <= g.lon[-1]))
                t.insert('end', f'  puntos fuera de su cobertura: {fuera}\n')
        t.insert('end', '\nLos valores de N del mapa son los del archivo, sin corrección de marea. '
                        'Coordenadas geodésicas (grados); el mapa conserva la proporción real '
                        '(1° de longitud = cos φ · 1° de latitud).\n')
        for x in info:
            t.insert('end', x + '\n')
        t.configure(state='disabled')
        self.estado.set('Mapa actualizado.' if not self.mapa_cargando else 'Leyendo mallas para el mapa…')

    def _mapa_cursor(self, ev):
        if not self.mapa_geo:
            return
        ox, oy, pw, ph, lon0, lon1, lat0, lat1, mf = self.mapa_geo
        if not (ox <= ev.x <= ox + pw and oy <= ev.y <= oy + ph):
            return
        lo = lon0 + (ev.x - ox) / pw * (lon1 - lon0)
        la = lat1 - (ev.y - oy) / ph * (lat1 - lat0)
        txt = f'lat {la:8.4f}°   lon {lo:9.4f}°'
        if mf is not None:
            try:
                z, _, _ = mp.interpolar(mf, 'Bilineal', la, lo, None, None, None)
                txt += f'   N {z:8.3f} m'
            except mp.NoAplicable:
                txt += '   N —'
        self.lbl_cursor.configure(text=txt)

    # ---------------------------------------------------------- 5. Calcular
    def _tab_calcular(self):
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  5 · Calcular  ')
        self._paso(fr, 'Paso 3 — Parámetros y ejecución',
                   'Los valores por defecto son los de appQGeoidCOL2023_0.m. Presione «Ejecutar cálculo»; '
                   'con dos mallas tarda alrededor de un minuto. La bitácora muestra cada paso.')
        top = ttk.Frame(fr)
        top.pack(fill='x', padx=14)
        par = ttk.LabelFrame(top, text='Parámetros')
        par.pack(side='left', fill='y')
        self.v = {'idw_p': tk.StringVar(value='2'), 'idw_k': tk.StringVar(value='16'),
                  'krig_n': tk.StringVar(value='16'), 'krig_ventana': tk.StringVar(value='6'),
                  'umbral_mad': tk.StringVar(value='5'), 'modelo': tk.StringVar(value='18'),
                  'salida': tk.StringVar(value=os.path.join(AQUI, 'resultados'))}
        filas = [('IDW: exponente p', 'idw_p', '0.5 – 6; mayor p = más peso a los nodos cercanos.  '
                                               '«Sugerir» busca el mejor p'),
                 ('IDW: número de nodos k', 'idw_k', '3 – 64'),
                 ('Kriging: número de vecinos n', 'krig_n', '4 – 64'),
                 ('Kriging: semiancho de ventana', 'krig_ventana', '3 – 12 nodos (para ajustar el variograma)'),
                 ('Depuración: umbral de error grueso', 'umbral_mad', 'veces la dispersión típica; mayor = aparta menos puntos')]
        for r, (t, k, ay) in enumerate(filas):
            ttk.Label(par, text=t).grid(row=r, column=0, sticky='w', padx=8, pady=3)
            ce = ttk.Frame(par)
            ce.grid(row=r, column=1, sticky='w')
            ttk.Entry(ce, textvariable=self.v[k], width=7).pack(side='left')
            if k == 'idw_p':
                ttk.Button(ce, text='Sugerir', style='Mini.TButton',
                           command=self._sugerir_p).pack(side='left', padx=(4, 0))
            ttk.Label(par, text=ay, style='Ayuda.TLabel').grid(row=r, column=2, sticky='w', padx=8)
        ttk.Label(par, text='Modelo de ajuste').grid(row=5, column=0, sticky='nw', padx=8, pady=3)
        fm = ttk.Frame(par)
        fm.grid(row=5, column=1, columnspan=2, sticky='w')
        ttk.Radiobutton(fm, text=NOMBRE_MODELO['18'] + '  (δH, x₁, x₂ — recomendado)', value='18',
                        variable=self.v['modelo']).pack(anchor='w')
        ttk.Radiobutton(fm, text=NOMBRE_MODELO['23'] + '  (x₀ … x₃)', value='23',
                        variable=self.v['modelo']).pack(anchor='w')
        ttk.Label(par, text='Carpeta de resultados').grid(row=6, column=0, sticky='w', padx=8, pady=(8, 6))
        fs = ttk.Frame(par)
        fs.grid(row=6, column=1, columnspan=2, sticky='we', pady=(8, 6))
        ttk.Entry(fs, textvariable=self.v['salida'], width=52).pack(side='left')
        ttk.Button(fs, text='…', width=3, command=self._elegir_salida).pack(side='left', padx=4)

        acc = ttk.Frame(top)
        acc.pack(side='left', fill='both', expand=True, padx=(16, 0))
        self.btn_ejecutar = ttk.Button(acc, text='▶  Ejecutar cálculo', style='Accion.TButton',
                                       command=self._ejecutar)
        self.btn_ejecutar.pack(anchor='w', pady=(14, 6))
        self.btn_detener = ttk.Button(acc, text='Detener', command=self._detener, state='disabled')
        self.btn_detener.pack(anchor='w')
        self.pb = ttk.Progressbar(acc, maximum=1.0, length=420, style='green.Horizontal.TProgressbar')
        self.pb.pack(anchor='w', pady=(14, 2))
        self.lbl_prog = ttk.Label(acc, text='', style='Ayuda.TLabel')
        self.lbl_prog.pack(anchor='w')

        bit = ttk.LabelFrame(fr, text='Bitácora')
        bit.pack(fill='both', expand=True, padx=14, pady=10)
        self.txt_log = tk.Text(bit, wrap='none', font=('Consolas', 10), bg='#1b2631', fg='#e8eef4',
                               insertbackground=BLANCO, relief='flat', padx=8, pady=6)
        sy = ttk.Scrollbar(bit, command=self.txt_log.yview)
        sx = ttk.Scrollbar(bit, orient='horizontal', command=self.txt_log.xview)
        self.txt_log.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.pack(side='right', fill='y')
        sx.pack(side='bottom', fill='x')
        self.txt_log.pack(fill='both', expand=True)
        self.txt_log.tag_configure('mejor', foreground='#58d68d')
        self.txt_log.tag_configure('aviso', foreground='#f5b041')

    def _sugerir_p(self):
        """Ventana para sugerir el exponente p de IDW por validación cruzada."""
        if not self.mallas:
            messagebox.showinfo('Sugerir p', 'Primero agregue al menos una malla (pestaña 3).')
            return
        try:
            k = self._leer_param('idw_k', int, 3, 64, 'IDW k')
        except ValueError as e:
            messagebox.showwarning('Sugerir p', str(e))
            return
        win = tk.Toplevel(self.root)
        win.title('Sugerir el exponente p de IDW')
        win.configure(bg=FONDO)
        win.transient(self.root)
        ttk.Label(win, text='¿Qué exponente p interpola mejor la malla?', style='Paso.TLabel'
                  ).pack(anchor='w', padx=14, pady=(12, 2))
        ttk.Label(win, text='Se retira la mitad de los nodos de la malla y, con IDW, se predicen 400 de los nodos '
                            'retirados para cada p entre 0.5 y 6 (con k = ' + str(k) + ' nodos). Se compara lo '
                            'predicho con el valor real del archivo; el p con menor error (RMSE) es el sugerido. '
                            'Es la misma validación cruzada de appQGeoidCOL2023_0.m.',
                  style='Ayuda.TLabel', wraplength=640, justify='left').pack(anchor='w', padx=14)
        fila = ttk.Frame(win)
        fila.pack(fill='x', padx=14, pady=8)
        ttk.Label(fila, text='Malla:').pack(side='left')
        cb = ttk.Combobox(fila, state='readonly', width=26, values=[m['nombre'] for m in self.mallas])
        cb.current(0)
        cb.pack(side='left', padx=6)
        btn = ttk.Button(fila, text='Calcular', style='Mini.TButton')
        btn.pack(side='left', padx=6)
        cv = tk.Canvas(win, width=640, height=280, bg=BLANCO, highlightthickness=1,
                       highlightbackground='#d5dbe3')
        cv.pack(padx=14)
        lbl = ttk.Label(win, text='Elija la malla y presione «Calcular».', style='Ayuda.TLabel')
        lbl.pack(anchor='w', padx=14, pady=6)
        usar = ttk.Button(win, text='Usar el p sugerido', style='Accion.TButton', state='disabled')
        usar.pack(anchor='e', padx=14, pady=(0, 12))
        estado = {}

        def terminado(nombre, res):
            if not win.winfo_exists():
                return
            btn.configure(state='normal')
            if isinstance(res, str):
                lbl.configure(text=res)
                return
            p_opt, tabla, n = res
            estado['p'] = p_opt
            rm = {p: r for p, r, _ in tabla}
            texto = (f'{nombre}: p sugerido = {p_opt:g}  (RMSE = {rm[p_opt]:.1f} mm, {n} nodos de '
                     f'control). Con p = 2 el RMSE es {rm.get(2.0, float("nan")):.1f} mm.')
            if p_opt >= max(rm):
                texto += (' El error sigue bajando hasta el máximo permitido (p = 6): en esta malla conviene '
                          'dar mucho más peso a los nodos más cercanos.')
            lbl.configure(text=texto)
            usar.configure(state='normal', text=f'Usar p = {p_opt:g}')
            self._grafico_p(cv, tabla, p_opt)

        def calcular():
            m = self.mallas[cb.current()]
            btn.configure(state='disabled')
            usar.configure(state='disabled')
            lbl.configure(text='Calculando… (con mallas grandes tarda unos segundos)')
            cv.delete('all')

            def trabajo():
                try:
                    g = self.mapa_mallas.get(m['ruta']) or mp.cargar_malla(m['ruta'])
                    res = mp.validacion_idw(g, k)
                except SystemExit as e:
                    res = str(e)
                except Exception as e:
                    res = f'No se pudo calcular: {e}'
                self.q.put(('llamar', terminado, m['nombre'], res))
            threading.Thread(target=trabajo, daemon=True).start()

        def aplicar():
            if 'p' in estado:
                self.v['idw_p'].set(f'{estado["p"]:g}')
                self.estado.set(f'IDW: p = {estado["p"]:g} (sugerido por validación cruzada).')
                win.destroy()

        btn.configure(command=calcular)
        usar.configure(command=aplicar)

    def _grafico_p(self, cv, tabla, p_opt):
        cv.delete('all')
        W, H = int(cv.cget('width')), int(cv.cget('height'))
        L, R, T, B = 70, 20, 20, 46
        ps = [t[0] for t in tabla]
        rs = [t[1] for t in tabla]
        y0, y1 = min(rs), max(rs)
        m = (y1 - y0) * 0.1 or 1.0
        y0, y1 = y0 - m, y1 + m
        X = lambda v: L + (v - ps[0]) / (ps[-1] - ps[0]) * (W - L - R)
        Y = lambda v: H - B - (v - y0) / (y1 - y0) * (H - T - B)
        for v in numeros_bonitos(y0, y1, 5):
            if y0 <= v <= y1:
                cv.create_line(L, Y(v), W - R, Y(v), fill='#eef1f5')
                cv.create_text(L - 6, Y(v), text=fmt_num(v), anchor='e', font=('Segoe UI', 9))
        for v in range(1, 7):
            cv.create_text(X(v), H - B + 14, text=str(v), font=('Segoe UI', 9))
        cv.create_rectangle(L, T, W - R, H - B, outline='#9aa5b1')
        pts = []
        for p, r in zip(ps, rs):
            pts += [X(p), Y(r)]
        cv.create_line(*pts, fill=AZUL2, width=2)
        for p, r in zip(ps, rs):
            cv.create_oval(X(p) - 3, Y(r) - 3, X(p) + 3, Y(r) + 3, fill=AZUL2, outline='')
        r_opt = dict(zip(ps, rs))[p_opt]
        cv.create_line(X(p_opt), T, X(p_opt), H - B, fill=VERDE, dash=(4, 3))
        cv.create_oval(X(p_opt) - 6, Y(r_opt) - 6, X(p_opt) + 6, Y(r_opt) + 6, fill=VERDE, outline=BLANCO)
        cv.create_text(X(p_opt) + 8, T + 10, text=f'p sugerido = {p_opt:g}', anchor='w', fill=VERDE,
                       font=('Segoe UI', 10, 'bold'))
        cv.create_text((L + W - R) / 2, H - 14, text='exponente p', font=('Segoe UI', 10))
        cv.create_text(16, (T + H - B) / 2, text='RMSE [mm]', angle=90, font=('Segoe UI', 10))

    def _elegir_salida(self):
        r = filedialog.askdirectory(title='Carpeta de resultados', initialdir=self.v['salida'].get())
        if r:
            self.v['salida'].set(r)

    def _leer_param(self, k, tipo, lo, hi, nombre):
        try:
            x = tipo(self.v[k].get().replace(',', '.'))
        except ValueError:
            raise ValueError(f'{nombre}: «{self.v[k].get()}» no es un número válido.')
        if not lo <= x <= hi:
            raise ValueError(f'{nombre} debe estar entre {lo} y {hi}.')
        return x

    def _armar_config(self):
        """Valida lo que el usuario eligió en los pasos 1-3; devuelve cfg o lanza ValueError."""
        ruta = self.v_puntos.get().strip()
        if not ruta or not os.path.isfile(ruta):
            raise ValueError('Paso 1: elija un archivo de puntos válido.')
        if not self.encabezado or getattr(self, 'ruta_leida', None) != os.path.abspath(ruta):
            self._cargar_puntos(ruta)
            raise ValueError('Paso 1: se leyó el archivo de puntos indicado. Revise las columnas '
                             'sugeridas y vuelva a presionar «Ejecutar cálculo».')
        idx = {}
        for k, cb in self.cb_col.items():
            idx[k] = None if cb.get() in ('', NO_VIENE) else self.encabezado.index(cb.get())
        if idx['lat'] is None or idx['lon'] is None:
            raise ValueError('Paso 1: indique las columnas de latitud y longitud.')
        if idx['h'] is None and idx['H'] is None:
            raise ValueError('Paso 1: indique al menos una altura (h o H).')
        usadas = [v for v in idx.values() if v is not None]
        if len(usadas) != len(set(usadas)):
            raise ValueError('Paso 1: hay una misma columna asignada a dos datos distintos.')
        sis_h = self._sis(self.cb_sis_h) if idx['h'] is not None else None
        sis_H = self._sis(self.cb_sis_H) if idx['H'] is not None else None
        if idx['h'] is not None and sis_h is None:
            raise ValueError('Paso 1: elija el sistema de mareas de h.')
        if idx['H'] is not None and sis_H is None:
            raise ValueError('Paso 1: elija el sistema de mareas de H.')
        mallas = []
        if idx['h'] is not None and idx['H'] is not None:
            for m in self.mallas:
                if not m['sistema']:
                    raise ValueError(f'Paso 2: elija el sistema de mareas de la malla {m["nombre"]}.')
                if not os.path.isfile(m['ruta']):
                    raise ValueError(f'Paso 2: no existe la malla {m["ruta"]}.')
                mallas.append((m['ruta'], None, m['sistema']))
        cfg = types.SimpleNamespace(
            puntos=ruta, idx=idx, sis_h=sis_h, sis_H=sis_H, mallas=mallas,
            idw_p=self._leer_param('idw_p', float, 0.5, 6, 'IDW p'),
            idw_k=self._leer_param('idw_k', int, 3, 64, 'IDW k'),
            krig_n=self._leer_param('krig_n', int, 4, 64, 'Kriging n'),
            krig_ventana=self._leer_param('krig_ventana', int, 3, 12, 'Kriging semiancho'),
            umbral_mad=self._leer_param('umbral_mad', float, 1, 100, 'Umbral R6'),
            modelo=self.v['modelo'].get(), salida=self.v['salida'].get().strip() or
            os.path.join(AQUI, 'resultados'))
        return cfg

    def _ejecutar(self):
        if self.hilo and self.hilo.is_alive():
            return
        try:
            cfg = self._armar_config()
        except ValueError as e:
            messagebox.showwarning('Falta información', str(e))
            return
        if not cfg.mallas and cfg.idx['h'] is not None and cfg.idx['H'] is not None:
            if not messagebox.askyesno('Sin mallas', 'No agregó ninguna malla de N (paso 2). Solo se '
                                       'depurarán y exportarán los puntos en zero tide.\n\n¿Continuar?'):
                return
        self.txt_log.delete('1.0', 'end')
        self.pb['value'] = 0
        self.detener = False
        self.btn_ejecutar.configure(state='disabled')
        self.btn_detener.configure(state='normal')
        self.estado.set('Calculando…')
        mp.LOG_HOOK = lambda t: self.q.put(('log', t))

        def trabajo():
            try:
                r = mp.ejecutar(cfg, lambda f, t: self.q.put(('prog', f, t)), lambda: self.detener)
                self.q.put(('fin', r))
            except mp.Cancelado:
                self.q.put(('cancelado', None))
            except SystemExit as e:
                self.q.put(('error', str(e)))
            except Exception as e:      # se muestra al usuario, no se oculta
                import traceback
                self.q.put(('error', f'{e}\n\n{traceback.format_exc()}'))

        self.hilo = threading.Thread(target=trabajo, daemon=True)
        self.hilo.start()

    def _detener(self):
        self.detener = True
        self.estado.set('Deteniendo…')

    def _leer_cola(self):
        try:
            while True:
                msg = self.q.get_nowait()
                tipo = msg[0]
                if tipo == 'log':
                    t = msg[1]
                    tag = 'mejor' if ('MEJOR' in t) else 'aviso' if ('AVISO' in t or 'ADVERTENCIA' in t) else ''
                    self.txt_log.insert('end', t + '\n', tag)
                    self.txt_log.see('end')
                elif tipo == 'prog':
                    self.pb['value'] = msg[1]
                    self.lbl_prog.configure(text=msg[2])
                elif tipo == 'llamar':
                    msg[1](*msg[2:])
                elif tipo == 'mapa_ok':
                    self.mapa_mallas[msg[1]] = msg[2]
                    self.mapa_cargando.discard(msg[1])
                    self._dibujar_mapa()
                elif tipo == 'mapa_error':
                    self.mapa_cargando.discard(msg[1])
                    messagebox.showerror('Mapa', msg[2])
                else:
                    self.btn_ejecutar.configure(state='normal')
                    self.btn_detener.configure(state='disabled')
                    mp.LOG_HOOK = None
                    if tipo == 'fin':
                        self.salida = msg[1]
                        self._mostrar_resultados()
                        self._preparar_kriging()
                        self._dibujar_mapa()
                        self.estado.set(f'Terminado. Resultados en {self.salida["carpeta"]}')
                        if self.salida['resultados']:
                            self.nb.select(5)
                        elif self.salida['aviso']:
                            messagebox.showinfo('Resultado', self.salida['aviso'])
                    elif tipo == 'cancelado':
                        self.estado.set('Cálculo detenido por el usuario.')
                    else:
                        self.estado.set('Error: revise el mensaje.')
                        messagebox.showerror('Error en el cálculo', msg[1][:3000])
        except queue.Empty:
            pass
        self.root.after(100, self._leer_cola)

    # -------------------------------------------------------- 5. Resultados
    def _tab_resultados(self):
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  6 · Resultados  ')
        self._paso(fr, 'Resultados — modelo de ajuste y mejor método',
                   'Para cada malla: los parámetros del modelo de ajuste con su desviación estándar, σ₀, R² y R² '
                   'ajustado de cada método de interpolación. En verde, el método con mayor R² ajustado.')
        top = ttk.Frame(fr)
        top.pack(fill='x', padx=14)
        self.lbl_dep = ttk.Label(top, text='Aún no hay resultados: ejecute el cálculo en la pestaña 5.')
        self.lbl_dep.pack(side='left')
        ttk.Button(top, text='Abrir carpeta de resultados', command=self._abrir_carpeta).pack(side='right')
        self.nb_res = ttk.Notebook(fr)
        self.nb_res.pack(fill='both', expand=True, padx=14, pady=8)

    def _abrir_carpeta(self):
        carpeta = self.salida['carpeta'] if self.salida else self.v['salida'].get()
        if not os.path.isdir(carpeta):
            messagebox.showinfo('Resultados', 'La carpeta todavía no existe.')
            return
        try:
            if sys.platform.startswith('win'):
                os.startfile(carpeta)                       # noqa: B606 (solo Windows)
            elif sys.platform == 'darwin':
                subprocess.Popen(['open', carpeta])
            else:
                subprocess.Popen(['xdg-open', carpeta])
        except OSError as e:
            messagebox.showinfo('Resultados', f'Carpeta: {carpeta}\n({e})')

    def _mostrar_resultados(self):
        s = self.salida
        for t in self.nb_res.tabs():
            self.nb_res.forget(t)
        n_rech = len(s['rechazos'])
        self.lbl_dep.configure(text=f'Depuración: {len(s["puntos"])} puntos aceptados, {n_rech} rechazados '
                                    f'(ver pestaña «Puntos rechazados»).')
        for r in s['resultados']:
            self._pestana_malla(r)
        self._pestana_rechazados(s)
        f = ttk.Frame(self.nb_res)
        self.nb_res.add(f, text='Archivos generados')
        t = tk.Text(f, height=12, font=('Consolas', 10), bg=BLANCO, relief='flat', padx=10, pady=8)
        t.pack(fill='both', expand=True)
        t.insert('end', 'Carpeta: ' + s['carpeta'] + '\n\n')
        expl = {'coordenadas_zero_tide.xlsx': 'Excel: puntos aceptados con h y H originales, corrección y zero tide',
                'coordenadas_zero_tide.csv': 'lo mismo en texto CSV',
                'puntos_rechazados.xlsx': 'Excel: puntos apartados en la depuración y el motivo',
                'puntos_rechazados.csv': 'lo mismo en texto CSV',
                'resumen_modelos.xlsx': 'Excel: parámetros, σ₀, R² y R² ajustado de todos los métodos y mallas',
                'resumen_modelos.csv': 'lo mismo en texto CSV',
                'reporte.txt': 'bitácora completa'}
        for nombre in sorted(os.listdir(s['carpeta'])):
            e = expl.get(nombre, '')
            if nombre.startswith('interpolacion_'):
                e = 'N, l = h−H−N y residuos v de los 6 métodos en cada punto; σ y modelo del kriging'
            elif nombre.startswith('modelo_parametrico_'):
                e = 'parámetros del modelo de ajuste para esa malla'
            elif nombre.endswith('_zero_tide.txt'):
                e = 'malla convertida a zero tide (mismo formato del original)'
            t.insert('end', f'  {nombre:<42} {e}\n')
        t.configure(state='disabled')
        self.nb_res.select(0)

    def _pestana_malla(self, r):
        m = r['malla']
        f = ttk.Frame(self.nb_res)
        self.nb_res.add(f, text=f'{m.nombre}  ({SIS_TEXTO.get(m.sistema)})')
        a_mejor = r['ajustes'][r['mejor']]
        txt = (f'Mejor método: {r["mejor"]}    R²_adj = {a_mejor["R2adj"]:.6f}    '
               f'δH = {a_mejor["x"][0]:.4f} ± {a_mejor["sx"][0]:.4f} m    σ₀ = {a_mejor["sigma0"]:.4f} m    '
               f'n = {a_mejor["n"]} puntos')
        ttk.Label(f, text=txt, style='Mejor.TLabel').pack(fill='x', pady=(8, 4))
        if r['menor_s0'] != r['mejor']:
            ttk.Label(f, text=f'Nota: con el criterio de menor σ₀ el mejor sería {r["menor_s0"]}. '
                              f'Las diferencias entre métodos suelen ser pequeñas; compare ambos criterios.',
                      style='Ayuda.TLabel').pack(anchor='w')
        cuerpo = ttk.Frame(f)
        cuerpo.pack(fill='both', expand=True, pady=6)
        nx = r['nombres_x']
        etiquetas = {'dH': 'δH [m]', 'x0': 'x₀ [m]', 'x1': 'x₁', 'x2': 'x₂', 'x3': 'x₃'}
        cols = ['metodo'] + nx + ['s0', 'r2', 'r2a']
        tv = ttk.Treeview(cuerpo, columns=cols, show='headings', height=7)
        for c, t in zip(cols, ['Método'] + [etiquetas[x] for x in nx] + ['σ₀ [m]', 'R²', 'R² ajustado']):
            tv.heading(c, text=t)
            tv.column(c, width=170 if c == 'metodo' else 135, anchor='w' if c == 'metodo' else 'e')
        tv.tag_configure('mejor', background=VERDE_CLARO, foreground=VERDE)
        for met in mp.METODOS:
            a = r['ajustes'][met]
            vals = [met] + [f'{x:.4f} ± {sx:.4f}' for x, sx in zip(a['x'], a['sx'])] + \
                   [f'{a["sigma0"]:.4f}', f'{a["R2"]:.6f}', f'{a["R2adj"]:.6f}']
            tv.insert('', 'end', values=vals, tags=('mejor',) if met == r['mejor'] else ())
        tv.pack(side='left', fill='both', expand=True)
        cv = tk.Canvas(cuerpo, width=380, height=230, bg=BLANCO, highlightthickness=0)
        cv.pack(side='left', fill='y', padx=(10, 0))
        cv.bind('<Configure>', lambda e, cv=cv, r=r: self._barras_r2(cv, r))
        unidades = ('x₁ en m/° (norte–sur), x₂ en m/° (este–oeste). φ₀ = {:.4f}°, λ₀ = {:.4f}°.'
                    if r['nombres_x'][0] == 'dH' else 'x₀…x₃ en m. φ₀ = {:.4f}°, λ₀ = {:.4f}°.')
        ttk.Label(f, text=unidades.format(a_mejor['phi0'], a_mejor['lam0']) +
                  '  Ajuste por mínimos cuadrados con todos los puntos con el mismo peso.',
                  style='Ayuda.TLabel').pack(anchor='w')

    def _barras_r2(self, cv, r):
        cv.delete('all')
        W, H = cv.winfo_width(), cv.winfo_height()
        if W < 50 or H < 50:
            return
        vals = [r['ajustes'][m]['R2adj'] for m in mp.METODOS]
        vmax = max(max(vals), 1e-9)
        cv.create_text(10, 12, text='R² ajustado por método', anchor='w', fill=AZUL,
                       font=('Segoe UI', 10, 'bold'))
        x0, alto = 132, (H - 36) / len(vals)
        for k, (met, v) in enumerate(zip(mp.METODOS, vals)):
            y = 30 + k * alto
            largo = max(v, 0) / vmax * (W - x0 - 70)
            color = VERDE if met == r['mejor'] else AZUL2
            cv.create_text(x0 - 6, y + alto / 2, text=met, anchor='e', font=('Segoe UI', 9))
            cv.create_rectangle(x0, y + 3, x0 + largo, y + alto - 3, fill=color, outline='')
            cv.create_text(x0 + largo + 4, y + alto / 2, text=f'{v:.4f}', anchor='w',
                           font=('Segoe UI', 9, 'bold' if met == r['mejor'] else 'normal'))

    def _pestana_rechazados(self, s):
        f = ttk.Frame(self.nb_res)
        self.nb_res.add(f, text=f'Puntos rechazados ({len(s["rechazos"])})')
        resumen = '   '.join(f'{k}: {v}' for k, v in sorted(s['conteo_rechazos'].items()))
        ttk.Label(f, text=resumen or 'No se rechazó ningún punto.', style='Ayuda.TLabel',
                  wraplength=1150).pack(anchor='w', pady=6)
        cols = ('id', 'linea', 'lat', 'lon', 'h', 'H', 'regla', 'detalle')
        cont = ttk.Frame(f)
        cont.pack(fill='both', expand=True)
        tv = ttk.Treeview(cont, columns=cols, show='headings')
        for c, t, w in zip(cols, ('Id', 'Línea', 'Latitud', 'Longitud', 'h [m]', 'H [m]', 'Regla', 'Detalle'),
                           (60, 55, 95, 100, 85, 85, 210, 520)):
            tv.heading(c, text=t)
            tv.column(c, width=w, anchor='w')
        sb = ttk.Scrollbar(cont, command=tv.yview)
        tv.configure(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y')
        tv.pack(fill='both', expand=True)
        for r in s['rechazos']:
            tv.insert('', 'end', values=[r.get(c) if r.get(c) is not None else '' for c in cols])

    # ----------------------------------------------- 6. Kriging/variograma
    def _tab_kriging(self):
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  7 · Kriging ordinario  ')
        self._paso(fr, 'Kriging ordinario — semivariograma y covariograma',
                   'El semivariograma se ajusta en una ventana local (uno por celda de la malla que contiene '
                   'puntos). No hace falta revisarlos uno por uno: la vista «Resumen» los reúne todos en una sola '
                   'gráfica. Con «Un punto» se ve el de un punto: puntos negros = valores experimentales (más '
                   'grandes = más pares de nodos), líneas = las 4 curvas ajustadas, línea gruesa = la elegida '
                   '(menor SSRw).')
        ctrl = ttk.Frame(fr)
        ctrl.pack(fill='x', padx=14)
        self.v_vista = tk.StringVar(value='resumen')
        ttk.Label(ctrl, text='Vista:').pack(side='left')
        ttk.Radiobutton(ctrl, text='Resumen de todos', value='resumen', variable=self.v_vista,
                        command=self._dibujar_variograma).pack(side='left', padx=(4, 0))
        ttk.Radiobutton(ctrl, text='Un punto', value='punto', variable=self.v_vista,
                        command=self._dibujar_variograma).pack(side='left', padx=(6, 16))
        ttk.Label(ctrl, text='Malla:').pack(side='left')
        self.cb_k_malla = ttk.Combobox(ctrl, state='readonly', width=26)
        self.cb_k_malla.pack(side='left', padx=(4, 14))
        ttk.Label(ctrl, text='Punto (id):').pack(side='left')
        self.cb_k_punto = ttk.Combobox(ctrl, width=14)
        self.cb_k_punto.pack(side='left', padx=(4, 14))
        self.v_tipo = tk.StringVar(value='semi')
        ttk.Radiobutton(ctrl, text='Semivariograma γ(h)', value='semi', variable=self.v_tipo,
                        command=self._dibujar_variograma).pack(side='left')
        ttk.Radiobutton(ctrl, text='Covariograma C(h)', value='cov', variable=self.v_tipo,
                        command=self._dibujar_variograma).pack(side='left', padx=(6, 14))
        self.v_todos = tk.BooleanVar(value=True)
        ttk.Checkbutton(ctrl, text='Mostrar los 4 modelos', variable=self.v_todos,
                        command=self._dibujar_variograma).pack(side='left')
        ttk.Button(ctrl, text='Guardar gráfico (.eps)', command=self._guardar_grafico).pack(side='right')
        self.cb_k_malla.bind('<<ComboboxSelected>>', lambda e: self._cambiar_malla_kriging())
        self.cb_k_punto.bind('<<ComboboxSelected>>', lambda e: self._ver_punto())
        self.cb_k_punto.bind('<Return>', lambda e: self._ver_punto())

        cuerpo = ttk.Frame(fr)
        cuerpo.pack(fill='both', expand=True, padx=14, pady=8)
        self.cv_var = tk.Canvas(cuerpo, bg=BLANCO, highlightthickness=1, highlightbackground='#d5dbe3')
        self.cv_var.pack(side='left', fill='both', expand=True)
        self.cv_var.bind('<Configure>', lambda e: self._dibujar_variograma())
        self.txt_var = tk.Text(cuerpo, width=46, font=('Consolas', 10), bg=BLANCO, relief='flat',
                               padx=10, pady=8, wrap='word')
        self.txt_var.pack(side='left', fill='y', padx=(10, 0))

    def _preparar_kriging(self):
        res = self.salida['resultados'] if self.salida else []
        self.cb_k_malla.configure(values=[r['malla'].nombre for r in res])
        if res:
            self.cb_k_malla.current(0)
        self._cambiar_malla_kriging()

    def _cambiar_malla_kriging(self):
        res = self.salida['resultados'] if self.salida else []
        r = next((x for x in res if x['malla'].nombre == self.cb_k_malla.get()), None)
        if r is not None:
            ids = [p['id'] for p in r['comunes']]
            self.cb_k_punto.configure(values=ids)
            if ids and self.cb_k_punto.get() not in ids:
                self.cb_k_punto.set(ids[0])
        self._dibujar_variograma()

    def _ver_punto(self):
        self.v_vista.set('punto')
        self._dibujar_variograma()

    def _datos_variograma(self):
        if not self.salida or not self.salida['resultados']:
            return None, 'Ejecute el cálculo (pestaña 5) para ver el semivariograma del kriging.'
        nombre = self.cb_k_malla.get()
        r = next((x for x in self.salida['resultados'] if x['malla'].nombre == nombre), None)
        if r is None:
            return None, 'Elija una malla.'
        pid = self.cb_k_punto.get().strip()
        info = r['sig_k'].get(pid)
        if info is None:
            return None, f'El punto «{pid}» no existe o fue rechazado en la depuración.'
        mod = r['variogramas'][info['celda']]
        p = next(q for q in r['comunes'] if q['id'] == pid)
        return (r, pid, info, mod, p), None

    def _dibujar_variograma(self):
        if self.v_vista.get() == 'resumen' and self.salida and self.salida['resultados']:
            self._dibujar_resumen()
        else:
            self._dibujar_punto()

    def _dibujar_resumen(self):
        """Todos los semivariogramas (o covariogramas) locales en una gráfica: mediana y franja P10-P90."""
        cv = self.cv_var
        cv.delete('all')
        W, H = cv.winfo_width(), cv.winfo_height()
        if W < 100 or H < 100:
            return
        r = next((x for x in self.salida['resultados'] if x['malla'].nombre == self.cb_k_malla.get()),
                 self.salida['resultados'][0])
        mods = list(r['variogramas'].values())
        if not mods:
            cv.create_text(W / 2, H / 2, text='No hay variogramas para esta malla.', fill=GRIS,
                           font=('Segoe UI', 12))
            return
        cov = self.v_tipo.get() == 'cov'
        fun = mp.modelo_covariograma if cov else mp.modelo_variograma

        def pct(v, q):
            v = sorted(v)
            return v[min(len(v) - 1, max(0, int(round(q * (len(v) - 1)))))]
        hmax = pct([m['lagMax'] for m in mods], 0.5)
        hs = [hmax * k / 80 for k in range(81)]
        curvas = [[fun(h if h > 0 else 1e-9, m['nugget'], m['sill'], m['rango'], m['tipo']) for h in hs]
                  for m in mods]
        p10 = [pct([c[k] for c in curvas], 0.10) for k in range(len(hs))]
        p50 = [pct([c[k] for c in curvas], 0.50) for k in range(len(hs))]
        p90 = [pct([c[k] for c in curvas], 0.90) for k in range(len(hs))]
        nb = 16                                       # mediana de los valores experimentales por rezago
        cajas = [[] for _ in range(nb)]
        for m in mods:
            for h, y in zip(m['hExp'], m['cExp'] if cov else m['gExp']):
                if h <= hmax:
                    cajas[min(nb - 1, int(h / hmax * nb))].append((h, y))
        exp = [(pct([c[0] for c in b], 0.5), pct([c[1] for c in b], 0.5)) for b in cajas if b]
        ys = p10 + p90 + [e[1] for e in exp] + [0.0]
        y0, y1 = min(ys), max(ys)
        mg = (y1 - y0) * 0.08 or 1e-6
        y0, y1 = (y0 - mg if y0 < 0 else 0.0), y1 + mg
        L, R, T, B = 90, 24, 46, 58
        X = lambda v: L + v / hmax * (W - L - R)
        Y = lambda v: H - B - (v - y0) / (y1 - y0) * (H - T - B)
        for v in numeros_bonitos(0, hmax):
            if v <= hmax * 1.0001:
                cv.create_line(X(v), T, X(v), H - B, fill='#eef1f5')
                cv.create_text(X(v), H - B + 14, text=fmt_num(v), font=('Segoe UI', 9), fill=NEGRO)
        for v in numeros_bonitos(y0, y1):
            if y0 <= v <= y1:
                cv.create_line(L, Y(v), W - R, Y(v), fill='#eef1f5')
                cv.create_text(L - 6, Y(v), text=fmt_num(v), anchor='e', font=('Segoe UI', 9), fill=NEGRO)
        banda = [c for h, y in zip(hs, p90) for c in (X(h), Y(y))] + \
                [c for h, y in zip(reversed(hs), reversed(p10)) for c in (X(h), Y(y))]
        cv.create_polygon(*banda, fill='#cfe0f3', outline='')
        for serie, estilo in ((p10, {'width': 1, 'dash': (4, 3)}), (p90, {'width': 1, 'dash': (4, 3)}),
                              (p50, {'width': 3})):
            pts = [c for h, y in zip(hs, serie) for c in (X(h), Y(y))]
            cv.create_line(*pts, fill=AZUL, **estilo)
        for h, y in exp:
            cv.create_oval(X(h) - 5, Y(y) - 5, X(h) + 5, Y(y) + 5, fill=NEGRO, outline=BLANCO)
        cv.create_rectangle(L, T, W - R, H - B, outline='#9aa5b1')
        if y0 < 0:
            cv.create_line(L, Y(0), W - R, Y(0), fill='#9aa5b1', dash=(2, 2))
        cv.create_text((L + W - R) / 2, H - 18, text='distancia h  [km]', font=('Segoe UI', 10), fill=NEGRO)
        cv.create_text(18, (T + H - B) / 2, text='C(h)  [m²]' if cov else 'γ(h)  [m²]', angle=90,
                       font=('Segoe UI', 10), fill=NEGRO)
        cv.create_text(L, 22, anchor='w', fill=AZUL, font=('Segoe UI', 12, 'bold'),
                       text=f'{"Covariogramas" if cov else "Semivariogramas"} locales — {r["malla"].nombre}: '
                            f'resumen de {len(mods)}')
        lx, ly = (W - R - 300, T + 14) if cov else (L + 14, T + 14)
        cv.create_rectangle(lx - 8, ly - 12, lx + 290, ly + 52, fill=BLANCO, outline='#d5dbe3')
        cv.create_line(lx, ly, lx + 24, ly, fill=AZUL, width=3)
        cv.create_text(lx + 32, ly, anchor='w', font=('Segoe UI', 9), text='mediana de las curvas ajustadas')
        cv.create_rectangle(lx, ly + 14, lx + 24, ly + 26, fill='#cfe0f3', outline='')
        cv.create_text(lx + 32, ly + 20, anchor='w', font=('Segoe UI', 9),
                       text='80 % central de las curvas (P10 – P90)')
        cv.create_oval(lx + 7, ly + 35, lx + 17, ly + 45, fill=NEGRO, outline='')
        cv.create_text(lx + 32, ly + 40, anchor='w', font=('Segoe UI', 9),
                       text='mediana de los valores experimentales')

        # panel de texto
        T_ = self.txt_var
        T_.configure(state='normal')
        T_.delete('1.0', 'end')
        T_.tag_configure('t', font=('Segoe UI', 11, 'bold'), foreground=AZUL)
        T_.tag_configure('e', foreground=VERDE, font=('Consolas', 10, 'bold'))
        n_pts = len(r['sig_k'])
        T_.insert('end', 'Resumen del kriging ordinario\n', 't')
        T_.insert('end', f'Semivariogramas ajustados: {len(mods)}\n'
                         f'Puntos interpolados:       {n_pts}\n'
                         'Los puntos de una misma celda comparten la ventana y, por tanto, el '
                         'semivariograma.\n\n')
        T_.insert('end', 'Modelo elegido (menor SSRw)\n', 't')
        cuenta = {}
        for m in mods:
            cuenta[m['tipo']] = cuenta.get(m['tipo'], 0) + 1
        for tipo, n in sorted(cuenta.items(), key=lambda t: -t[1]):
            T_.insert('end', f'{mp.NOMBRE_VARIOGRAMA[tipo]:<12} {n:>5}  ({100 * n / len(mods):.0f} %)\n',
                      'e' if n == max(cuenta.values()) else '')
        T_.insert('end', '\nValores típicos: mediana [P10 – P90]\n', 't')
        for nombre, clave, f in (('pepita C₀', 'nugget', '{:.2e} m²'), ('meseta parcial C', 'sill', '{:.2e} m²'),
                                 ('rango a', 'rango', '{:.1f} km')):
            v = [m[clave] for m in mods]
            T_.insert('end', f'{nombre:<17}{f.format(pct(v, .5))}\n{"":<17}[{f.format(pct(v, .1))} – '
                             f'{f.format(pct(v, .9))}]\n')
        tope = sum(1 for m in mods if m['rango'] >= 0.999 * m['lagMax'])
        T_.insert('end', f'\nRango en el máximo permitido: {100 * tope / len(mods):.0f} %\n')
        if tope:
            T_.insert('end', 'En esas ventanas N sigue cambiando con la distancia (tendencia regional), así '
                             'que la curva no llega a estabilizarse. No afecta mucho al kriging, que solo usa '
                             'los vecinos más cercanos.\n')
        sig = sorted(i['sigma'] * 1000 for i in r['sig_k'].values())
        T_.insert('end', '\nIncertidumbre del kriging (σ_KO)\n', 't')
        T_.insert('end', f'mediana {pct(sig, .5):.2f} mm\nP95     {pct(sig, .95):.2f} mm\n'
                         f'máximo  {sig[-1]:.2f} mm\n')
        if cov:
            T_.insert('end', '\nEn el covariograma, los valores experimentales se calculan respecto a la media '
                             'de cada ventana. Cuando N tiene tendencia regional quedan por debajo de la curva '
                             'del modelo: es esperable y no indica un error.\n')
        T_.insert('end', '\nPara ver el semivariograma de un punto, elíjalo en «Punto».\n')
        T_.configure(state='disabled')

    def _dibujar_punto(self):
        cv = self.cv_var
        cv.delete('all')
        W, H = cv.winfo_width(), cv.winfo_height()
        if W < 100 or H < 100:
            return
        datos, error = self._datos_variograma()
        self.txt_var.configure(state='normal')
        self.txt_var.delete('1.0', 'end')
        if error:
            cv.create_text(W / 2, H / 2, text=error, fill=GRIS, font=('Segoe UI', 11))
            self.txt_var.configure(state='disabled')
            return
        r, pid, info, mod, p = datos
        cov = self.v_tipo.get() == 'cov'
        hx = mod['hExp']
        yexp = mod['cExp'] if cov else mod['gExp']
        modelos = mod['todos'] if self.v_todos.get() else [t for t in mod['todos'] if t['tipo'] == mod['tipo']]
        fun = mp.modelo_covariograma if cov else mp.modelo_variograma
        hmax = mod['lagMax']
        muestras = [hmax * k / 160 for k in range(161)]
        curvas = []
        for t in modelos:
            ys = []
            for hh in muestras:
                hh2 = hh if hh > 0 else 1e-9          # γ(0+) = pepita (salto en el origen)
                ys.append(fun(hh2, t['nugget'], t['sill'], t['rango'], t['tipo']))
            curvas.append((t, ys))
        todos_y = list(yexp) + [y for _, ys in curvas for y in ys] + [0.0]
        ymin, ymax = min(todos_y), max(todos_y)
        margen = (ymax - ymin) * 0.08 or abs(ymax) * 0.1 or 1e-6
        ymin, ymax = (ymin - margen if ymin < 0 else 0.0), ymax + margen
        L, R, T, B = 86, 24, 46, 58
        X = lambda v: L + v / hmax * (W - L - R)
        Y = lambda v: H - B - (v - ymin) / (ymax - ymin) * (H - T - B)
        # rejilla y ejes
        for v in numeros_bonitos(0, hmax):
            if v <= hmax * 1.0001:
                cv.create_line(X(v), T, X(v), H - B, fill='#eef1f5')
                cv.create_text(X(v), H - B + 14, text=fmt_num(v), font=('Segoe UI', 9), fill=GRIS)
        for v in numeros_bonitos(ymin, ymax):
            if ymin <= v <= ymax:
                cv.create_line(L, Y(v), W - R, Y(v), fill='#eef1f5')
                cv.create_text(L - 6, Y(v), text=fmt_num(v), anchor='e', font=('Segoe UI', 9), fill=GRIS)
        cv.create_rectangle(L, T, W - R, H - B, outline='#9aa5b1')
        if ymin < 0:
            cv.create_line(L, Y(0), W - R, Y(0), fill='#9aa5b1', dash=(2, 2))
        cv.create_text((L + W - R) / 2, H - 18, text='distancia h  [km]  (proyección Origen-Nacional)',
                       font=('Segoe UI', 10), fill='#1b2631')
        cv.create_text(18, (T + H - B) / 2, text='C(h)  [m²]' if cov else 'γ(h)  [m²]', angle=90,
                       font=('Segoe UI', 10), fill='#1b2631')
        titulo = 'Covariograma' if cov else 'Semivariograma'
        cv.create_text(L, 22, anchor='w', fill=AZUL, font=('Segoe UI', 12, 'bold'),
                       text=f'{titulo} local — {r["malla"].nombre}, punto {pid}')
        # meseta (semivariograma) o C(0) (covariograma)
        elegido = next(t for t in mod['todos'] if t['tipo'] == mod['tipo'])
        meseta = elegido['nugget'] + elegido['sill']
        if ymin <= meseta <= ymax:
            cv.create_line(L, Y(meseta), W - R, Y(meseta), fill='#b03a2e', dash=(6, 4))
            cv.create_text(W - R - 4, Y(meseta) - 8 if not cov else Y(meseta) + 8, anchor='e',
                           fill='#b03a2e', font=('Segoe UI', 9),
                           text=('C(0) = C₀ + C' if cov else 'meseta C₀ + C') + f' = {fmt_num(meseta)} m²')
        # curvas
        for t, ys in curvas:
            es = t['tipo'] == mod['tipo']
            pts = []
            for hh, yy in zip(muestras, ys):
                pts += [X(hh), Y(min(max(yy, ymin), ymax))]
            estilo = {'width': 3} if es else {'width': 1.5, 'dash': (5, 3)}
            cv.create_line(*pts, fill=COLORES_MODELO[t['tipo']], **estilo)
        if cov:                                    # C(0) del modelo elegido (sin pepita: salto)
            cv.create_oval(X(0) - 4, Y(meseta) - 4, X(0) + 4, Y(meseta) + 4, fill='#b03a2e', outline='')
        # puntos experimentales
        npmax = max(mod['nPares'])
        for hh, yy, npar in zip(hx, yexp, mod['nPares']):
            rad = 3 + 5 * math.sqrt(npar / npmax)
            cv.create_oval(X(hh) - rad, Y(yy) - rad, X(hh) + rad, Y(yy) + rad, fill='#1b2631',
                           outline=BLANCO)
        # leyenda
        lx, ly = (L + 14, T + 14) if cov else (W - R - 240, H - B - 22 - 20 * len(curvas))
        cv.create_rectangle(lx - 8, ly - 12, lx + 230, ly + 20 * len(curvas) + 8, fill=BLANCO,
                            outline='#d5dbe3')
        cv.create_oval(lx, ly - 4, lx + 8, ly + 4, fill='#1b2631', outline='')
        cv.create_text(lx + 16, ly, anchor='w', font=('Segoe UI', 9), text='experimental (tamaño ∝ nº pares)')
        for k, (t, _) in enumerate(curvas):
            yy = ly + 20 * (k + 1)
            es = t['tipo'] == mod['tipo']
            estilo = {'width': 3} if es else {'width': 1.5, 'dash': (5, 3)}
            cv.create_line(lx, yy, lx + 24, yy, fill=COLORES_MODELO[t['tipo']], **estilo)
            cv.create_text(lx + 32, yy, anchor='w', font=('Segoe UI', 9, 'bold' if es else 'normal'),
                           text=f'{mp.NOMBRE_VARIOGRAMA[t["tipo"]]}' + ('  (elegido)' if es else ''))

        # panel de texto
        Nk = r['res']['Kriging ordinario'][pid]
        T_ = self.txt_var
        T_.tag_configure('t', font=('Segoe UI', 10, 'bold'), foreground=AZUL)
        T_.tag_configure('e', foreground=VERDE, font=('Consolas', 9, 'bold'))
        T_.insert('end', f'Punto {pid}\n', 't')
        T_.insert('end', f'φ = {p["lat"]:.6f}°   λ = {p["lon"]:.6f}°\n'
                         f'h_ZT = {p["h_zt"]:.3f} m   H_ZT = {p["H_zt"]:.3f} m\n\n')
        T_.insert('end', 'Resultado del kriging\n', 't')
        T_.insert('end', f'N (kriging)   = {Nk:.4f} m\n'
                         f'σ_KO          = {info["sigma"] * 1000:.2f} mm\n'
                         f'vecinos usados= {info["n_vecinos"]}\n'
                         f'nodos ventana = {mod["nNodos"]}\n'
                         f'l = h − H − N = {p["h_zt"] - p["H_zt"] - Nk:.4f} m\n\n')
        T_.insert('end', 'Modelo elegido (menor SSRw)\n', 't')
        T_.insert('end', f'{mp.NOMBRE_VARIOGRAMA[mod["tipo"]]}\n', 'e')
        T_.insert('end', f'pepita C₀       = {elegido["nugget"]:.3e} m²\n'
                         f'meseta parcial C= {elegido["sill"]:.3e} m²\n'
                         f'meseta C₀+C     = {meseta:.3e} m²\n'
                         f'rango a         = {elegido["rango"]:.2f} km\n'
                         f'varianza ventana= {mod["varianza"]:.3e} m²\n\n')
        T_.insert('end', 'Comparación de modelos (SSRw)\n', 't')
        for t in sorted(mod['todos'], key=lambda t: t['SSRw']):
            marca = '  ←' if t['tipo'] == mod['tipo'] else ''
            T_.insert('end', f'{mp.NOMBRE_VARIOGRAMA[t["tipo"]]:<12}{t["SSRw"]:.3e}  a={t["rango"]:.1f} km{marca}\n',
                      'e' if marca else '')
        T_.insert('end', '\nγ(h) = ½·prom[(N(x) − N(x+h))²]\nC(h) = C₀ + C − γ(h)\n'
                         'Ajuste: mín Σ nᵢ (γ̂ᵢ − γ(hᵢ))²\n')
        T_.configure(state='disabled')

    def _guardar_grafico(self):
        if not self.salida or not self.salida['resultados']:
            messagebox.showinfo('Gráfico', 'Primero ejecute el cálculo.')
            return
        tipo = 'covariograma' if self.v_tipo.get() == 'cov' else 'semivariograma'
        ini = f'{tipo}_{self.cb_k_malla.get()}_punto_{self.cb_k_punto.get()}.eps'
        r = filedialog.asksaveasfilename(title='Guardar gráfico', defaultextension='.eps',
                                         initialdir=self.salida['carpeta'], initialfile=ini,
                                         filetypes=[('PostScript encapsulado', '*.eps')])
        if r:
            self.cv_var.postscript(file=r, colormode='color')
            self.estado.set(f'Gráfico guardado: {r}')

    # ------------------------------------------------------------ inicio
    def _precargar(self):
        """Si los archivos del repositorio están junto al script, se cargan solos
        (los sistemas de mareas de los puntos siempre los elige el usuario)."""
        p = os.path.join(AQUI, 'coords_QGEOID2023.txt')
        if os.path.isfile(p):
            self.v_puntos.set(p)
            self._cargar_puntos(p)
        for f in ('QGeoidCOL2023.txt', 'coords_EGM2008.txt', 'EGM2008.txt'):
            ruta = os.path.join(AQUI, f)
            if os.path.isfile(ruta):
                self._agregar_malla(ruta)


def main():
    try:                                    # nitidez en pantallas de alta resolución (Windows)
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == '__main__':
    main()