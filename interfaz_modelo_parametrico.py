#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interfaz gráfica del modelo paramétrico GNSS/nivelación en Zero Tide.

Solo usa Python estándar (tkinter). Los cálculos los hace
modelo_parametrico_zero_tide.py (debe estar en la misma carpeta), así que la
interfaz y la línea de comandos dan exactamente los mismos resultados.

Ejecutar:   python interfaz_modelo_parametrico.py

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

import modelo_parametrico_zero_tide as mp

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
GRIS = '#5f6b7a'
VERDE = '#1e8449'
VERDE_CLARO = '#d5f5e3'
ROJO = '#c0392b'
COLORES_MODELO = {'sph': '#2e75b6', 'exp': '#e67e22', 'gau': '#8e44ad', 'cir': '#16a085'}

ELEGIR = '— elija —'
NO_VIENE = '(no viene)'
SIS_TEXTO = {'zero_tide': 'Zero tide', 'mean_tide': 'Mean tide', 'tide_free': 'Tide free'}
TEXTO_SIS = {v: k for k, v in SIS_TEXTO.items()}


def explicar_h(s):
    return {None: 'Elija el sistema de mareas de h.',
            'zero_tide': 'Ya está en zero tide: no se transforma.',
            'tide_free': 'Se aplica la Ec. (3) del paper:  h_ZT = h_TF + 0.0602 − 0.1790 sin²φ − 0.0018 sin⁴φ',
            'mean_tide': 'La corteza "mean tide" coincide con la "zero tide": no se transforma.'}[s]


def explicar_H(s):
    return {None: 'Elija el sistema de mareas de H.',
            'zero_tide': 'Ya está en zero tide: no se transforma.',
            'mean_tide': 'Se aplican las Ec. (1)-(2):  H_ZT = H_MT + ΔW/g,   ΔW/g ≈ −0.198 (3/2 sin²φ − 1/2)',
            'tide_free': 'H_ZT = H_TF + Δh(Ec. 3) − 0.3·[−0.198 (3/2 sin²φ − 1/2)]   (derivada, poco habitual)'}[s]


def explicar_N(s):
    return {None: 'Elija el sistema de mareas de la malla.',
            'zero_tide': 'Ya está en zero tide: no se transforma.',
            'tide_free': 'N_ZT = N_TF + 0.3·[−0.198 (3/2 sin²φ − 1/2)]   (derivada de la Ec. 2; Ekman 1989, k = 0.3). '
                         'Se exporta la malla convertida a <malla>_zero_tide.txt.',
            'mean_tide': 'N_ZT = N_MT + 0.198 (3/2 sin²φ − 1/2)   (derivada de la Ec. 2). '
                         'Se exporta la malla convertida a <malla>_zero_tide.txt.'}[s]


def fuente_de(nombre):
    n = nombre.lower()
    if 'egm' in n:
        return 'ICGEM (GFZ Potsdam)'
    if 'qgeoid' in n or 'geocol' in n:
        return 'IGAC – Red Geodésica'
    return '—'


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
        root.title('Modelo paramétrico GNSS/nivelación en Zero Tide — Guimarães et al. (2025)')
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
        self._tab_calcular()
        self._tab_resultados()
        self._tab_kriging()
        self.estado = tk.StringVar(value='Listo. Lea la guía (pestaña 1) y siga los pasos en orden.')
        tk.Label(root, textvariable=self.estado, anchor='w', bg=AZUL, fg=BLANCO,
                 font=('Segoe UI', 9), padx=10, pady=3).pack(fill='x', side='bottom')
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
        st.configure('Ayuda.TLabel', background=FONDO, foreground=GRIS, font=('Segoe UI', 9))
        st.configure('Formula.TLabel', background='#eaf1fb', foreground='#1b2631',
                     font=('Consolas', 10), padding=6)
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
        st.configure('Treeview', rowheight=24, font=('Segoe UI', 9))
        st.configure('Treeview.Heading', font=('Segoe UI', 9, 'bold'), foreground=AZUL)
        st.configure('green.Horizontal.TProgressbar', background=VERDE)

    def _encabezado(self):
        f = tk.Frame(self.root, bg=AZUL)
        f.pack(fill='x')
        tk.Label(f, text='Modelo paramétrico GNSS/nivelación en Zero Tide', bg=AZUL, fg=BLANCO,
                 font=('Segoe UI', 16, 'bold')).pack(anchor='w', padx=14, pady=(10, 0))
        tk.Label(f, text='Sistema de mareas → zero tide  ·  depuración  ·  6 métodos de interpolación  ·  '
                         'modelo de la Ec. (18)  ·  R² ajustado', bg=AZUL, fg='#d6e4f0',
                 font=('Segoe UI', 10)).pack(anchor='w', padx=14)
        enlace = tk.Label(f, text='Basado en: Guimarães et al. (2025), Earth, Planets and Space 77:98  '
                                  '(clic para abrir el paper)', bg=AZUL, fg='#aed6f1',
                          font=('Segoe UI', 9, 'underline'), cursor='hand2')
        enlace.pack(anchor='w', padx=14, pady=(0, 8))
        enlace.bind('<Button-1>', lambda e: webbrowser.open(URL_PAPER))

    def _paso(self, padre, titulo, ayuda):
        ttk.Label(padre, text=titulo, style='Paso.TLabel').pack(anchor='w', padx=14, pady=(12, 0))
        ttk.Label(padre, text=ayuda, style='Ayuda.TLabel', wraplength=1150,
                  justify='left').pack(anchor='w', padx=14, pady=(2, 8))

    # ------------------------------------------------------------- 1. Guía
    def _tab_guia(self):
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  1 · Guía y teoría  ')
        cont = ttk.Frame(fr)
        cont.pack(fill='both', expand=True, padx=12, pady=10)
        txt = tk.Text(cont, wrap='word', bg=BLANCO, relief='flat', padx=24, pady=18,
                      font=('Segoe UI', 10), spacing1=2, spacing3=4)
        sb = ttk.Scrollbar(cont, command=txt.yview)
        txt.configure(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y')
        txt.pack(side='left', fill='both', expand=True)
        txt.tag_configure('h1', font=('Segoe UI', 15, 'bold'), foreground=AZUL, spacing1=10, spacing3=6)
        txt.tag_configure('h2', font=('Segoe UI', 12, 'bold'), foreground=AZUL2, spacing1=12, spacing3=4)
        txt.tag_configure('f', font=('Consolas', 10), background='#eaf1fb', lmargin1=30, lmargin2=30)
        txt.tag_configure('nota', foreground=GRIS, font=('Segoe UI', 9, 'italic'))
        txt.tag_configure('b', font=('Segoe UI', 10, 'bold'))
        txt.tag_configure('lista', lmargin1=18, lmargin2=34)

        def link(texto, url):
            tag = f'link{len(txt.tag_names())}'
            txt.tag_configure(tag, foreground=AZUL2, underline=True)
            txt.tag_bind(tag, '<Button-1>', lambda e: webbrowser.open(url))
            txt.tag_bind(tag, '<Enter>', lambda e: txt.configure(cursor='hand2'))
            txt.tag_bind(tag, '<Leave>', lambda e: txt.configure(cursor=''))
            txt.insert('end', texto, tag)

        I = lambda t, *tags: txt.insert('end', t, tags)
        I('¿Qué hace este programa?\n', 'h1')
        I('Compara h − H con la ondulación del geoide N en puntos GNSS/nivelación y estima el '
          'desplazamiento del datum vertical local con el modelo paramétrico de ')
        link('Guimarães et al. (2025)', URL_PAPER)
        I('. Todo se lleva primero al mismo sistema de mareas (zero tide), como exige el paper.\n')
        I('Pasos (una pestaña por paso)\n', 'h2')
        for t in ('2 · Puntos: elija el archivo de puntos, diga qué columna es latitud, longitud, h y H, '
                  'y en qué sistema de mareas están h y H.',
                  '3 · Mallas: agregue las mallas de N (geoide/cuasigeoide) y diga su sistema de mareas.',
                  '4 · Calcular: revise los parámetros y presione «Ejecutar cálculo».',
                  '5 · Resultados: tabla de los 6 métodos, R² ajustado y mejor método; puntos rechazados.',
                  '6 · Kriging: semivariograma o covariograma del kriging en cualquier punto.'):
            I('  •  ' + t + '\n', 'lista')

        I('Puntos y mallas: no son lo mismo\n', 'h2')
        I('Puntos (coords_QGEOID2023.txt): ', 'b')
        I('son mediciones reales en campo, cada una con latitud φ, longitud λ, altura elipsoidal h '
          '(GNSS) y altura sobre el nivel del mar H (nivelación). Son los lugares donde se calcula el modelo.\n')
        I('Mallas (QGeoidCOL2023.txt, coords_EGM2008.txt): ', 'b')
        I('son modelos del geoide: una cuadrícula regular de nodos (cada 5′ o 2′) con el valor de N. '
          'Los puntos casi nunca caen sobre un nodo, por eso N se INTERPOLA desde los nodos de la malla '
          'hasta cada punto. Así, para cada punto se obtiene:\n')
        I('  l = h − H − N        (si todo fuera perfecto, l sería constante = desplazamiento del datum)\n', 'f')

        I('1. Sistema de mareas → zero tide (paper, sección 2.4)\n', 'h2')
        I('La marea permanente deforma la Tierra y el geoide. Hay tres convenciones: tide free, mean tide '
          'y zero tide. El paper trabaja en zero tide. Lo que ya está en zero tide NO se transforma.\n')
        I('  Ec. (1)  H_ZT = H_MT + ΔW_ZT / g\n'
          '  Ec. (2)  ΔW_ZT / g ≈ −0.198 (3/2 sin²φ − 1/2)                   [m]\n'
          '  Ec. (3)  h_ZT = h_TF + 0.0602 − 0.1790 sin²φ − 0.0018 sin⁴φ     [m]\n', 'f')
        I('Para las mallas de N (no está escrito en el paper; se deduce de la Ec. 2, Ekman 1989, k = 0.3):\n')
        I('  N_ZT = N_TF + 0.3 · [−0.198 (3/2 sin²φ − 1/2)]      N_ZT = N_MT + 0.198 (3/2 sin²φ − 1/2)\n', 'f')
        I('En Colombia, pasar EGM2008 de tide free a zero tide le suma a N entre 2.5 y 3.0 cm.\n', 'nota')

        I('Depuración de los puntos\n', 'h2')
        for t in ('R1 dato vacío o no numérico      R2 coordenada fuera de rango      R3 id o coordenadas duplicadas',
                  'R4 h = 0 con H > 1 m (altura elipsoidal faltante)      R5 |h − H| > 150 m (no plausible)',
                  'R6 error grosero: |l − mediana(l)| > umbral · σ robusta (σ = 1.4826·MAD; l con interpolación bilineal)',
                  'R7 punto fuera de la cobertura de la malla'):
            I('  •  ' + t + '\n', 'lista')
        I('Los rechazados se listan con su motivo en puntos_rechazados.csv.\n', 'nota')

        I('2. Interpolación de N: 6 métodos (como appQGeoidCOL2023_0.m)\n', 'h2')
        I('Todos calculan N = Σ wᵢ Nᵢ con Σ wᵢ = 1; ninguno extrapola fuera de la malla.\n')
        for t in ('Vecino más cercano: el nodo más próximo (w = 1).',
                  'Lineal por triángulos: la celda se parte por la diagonal SO–NE; pesos baricéntricos (3 nodos).',
                  'Bilineal: los 4 nodos de la celda; w = (1−u)(1−v), u(1−v), uv, (1−u)v.',
                  'Bicúbica: convolución cúbica de Keys (a = −0.5) con 4×4 = 16 nodos.',
                  'IDW: k nodos más cercanos, wᵢ ∝ 1/dᵢᵖ (distancias en Origen-Nacional, EPSG:9377).',
                  'Kriging ordinario: semivariograma local ajustado por mínimos cuadrados ponderados '
                  '(esférico, exponencial, gaussiano o circular) y sistema de kriging con n vecinos.'):
            I('  •  ' + t + '\n', 'lista')

        I('Semivariograma y covariograma (kriging)\n', 'h2')
        I('  γ(h) = ½ · promedio[ (N(x) − N(x+h))² ]          C(h) = C(0) − γ(h),   C(0) = C₀ + C\n'
          '  Esférico: γ = C₀ + C·(1.5 h/a − 0.5 (h/a)³)        Exponencial: γ = C₀ + C·(1 − e^(−h/a))\n'
          '  Gaussiano: γ = C₀ + C·(1 − e^(−(h/a)²))             C₀ = pepita, C = meseta parcial, a = rango\n', 'f')
        I('Se elige el modelo con menor suma de cuadrados ponderada por el número de pares de cada rezago. '
          'Luego se resuelve  [Γ 1; 1ᵀ 0]·[λ; μ] = [γ₀; 1]  y la varianza de kriging es σ²_KO = λᵀγ₀ + μ.\n')

        I('3. Modelo paramétrico, Ec. (18) del paper\n', 'h2')
        I('  l_P = h_P − H_P − N_P = δH + x₁ (φ_P − φ₀) + x₂ (λ_P − λ₀) cos φ_P\n'
          '  x̂ = (Aᵀ P A)⁻¹ Aᵀ P l,    P = I         (Ec. 19–20)\n', 'f')
        I('δH es el desplazamiento medio del datum vertical; x₁ y x₂ son las inclinaciones norte–sur y '
          'este–oeste; (φ₀, λ₀) es la posición media de los puntos. También está disponible la superficie '
          'de 4 parámetros de la Ec. (23):  x₀ + x₁ cosφ cosλ + x₂ cosφ sinλ + x₃ sinφ.\n')
        I('4. Coeficiente de determinación ajustado\n', 'h2')
        I('  R² = 1 − SSE/SST        R²_adj = 1 − (1 − R²)·(n − 1)/(n − u)       (u = nº de parámetros)\n', 'f')
        I('Se ajusta el modelo con el N de cada método; el método con MAYOR R²_adj es el mejor. Los 6 '
          'métodos se comparan con los mismos puntos. Como dato adicional se muestra σ₀ (el paper reporta '
          'desviaciones estándar).\n')

        I('Fuentes de los datos y referencia\n', 'h2')
        I('Paper:  ', 'b')
        I(CITA + '  ')
        link(URL_PAPER, URL_PAPER)
        I('\n')
        I('Malla EGM2008 (coords_EGM2008.txt):  ', 'b')
        I('ICGEM – GFZ Potsdam, ')
        link(URL_ICGEM, URL_ICGEM)
        I('\n')
        I('Malla QGeoidCOL2023 y puntos coords_QGEOID2023:  ', 'b')
        I('IGAC – Red Geodésica Nacional, ')
        link(URL_IGAC, URL_IGAC)
        I('\n')
        txt.configure(state='disabled')

    # ------------------------------------------------------------ 2. Puntos
    def _tab_puntos(self):
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  2 · Puntos  ')
        self._paso(fr, 'Paso 1 — Archivo de puntos GNSS/nivelación',
                   'Son los puntos medidos (ej. coords_QGEOID2023.txt del IGAC). En ellos se interpola N y se '
                   'calcula l = h − H − N. Formatos: texto con encabezado separado por «;», «,», tabulación '
                   'o espacios.')
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
                  ('lon', 'Longitud λ [grados]'), ('h', 'Altura elipsoidal h [m] (GNSS)'),
                  ('H', 'Altura sobre el nivel del mar H [m] (nivelación)')]
        for r, (k, t) in enumerate(textos):
            ttk.Label(col, text=t).grid(row=r, column=0, sticky='w', padx=8, pady=3)
            cb = ttk.Combobox(col, state='readonly', width=28)
            cb.grid(row=r, column=1, sticky='w', padx=8, pady=3)
            self.cb_col[k] = cb
        ttk.Label(col, text='Si una altura no viene en el archivo, elija «(no viene)». Para el modelo se '
                            'necesitan h y H.', style='Ayuda.TLabel', wraplength=420
                  ).grid(row=5, column=0, columnspan=2, sticky='w', padx=8, pady=(2, 6))

        mar = ttk.LabelFrame(abajo, text='¿En qué sistema de mareas vienen las alturas?')
        mar.pack(side='left', fill='both', expand=True)
        self.cb_sis_h = ttk.Combobox(mar, state='readonly', width=16, values=list(SIS_TEXTO.values()))
        self.cb_sis_H = ttk.Combobox(mar, state='readonly', width=16, values=list(SIS_TEXTO.values()))
        self.lbl_h = ttk.Label(mar, text=explicar_h(None), style='Formula.TLabel', wraplength=560)
        self.lbl_H = ttk.Label(mar, text=explicar_H(None), style='Formula.TLabel', wraplength=560)
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
                           lambda e: self.lbl_h.configure(text=explicar_h(self._sis(self.cb_sis_h))))
        self.cb_sis_H.bind('<<ComboboxSelected>>',
                           lambda e: self.lbl_H.configure(text=explicar_H(self._sis(self.cb_sis_H))))

    @staticmethod
    def _sis(cb):
        return TEXTO_SIS.get(cb.get())

    def _elegir_puntos(self):
        r = filedialog.askopenfilename(title='Archivo de puntos GNSS/nivelación', initialdir=AQUI,
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
        self._paso(fr, 'Paso 2 — Mallas de N (geoide / cuasigeoide)',
                   'Una malla es un MODELO del geoide: una cuadrícula regular de nodos (lat, lon, N). No son '
                   'puntos medidos. De cada malla se interpola N en los puntos del paso 1 con 6 métodos, y cada '
                   'malla se analiza por separado. Formatos: «lat lon N» sin encabezado (ej. QGeoidCOL2023.txt '
                   'del IGAC) o archivo ICGEM con encabezado (ej. coords_EGM2008.txt).')
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
        self.lbl_N = ttk.Label(sel, text=explicar_N(None), style='Formula.TLabel', wraplength=1100)
        self.lbl_N.grid(row=1, column=0, columnspan=2, sticky='we', padx=8, pady=(0, 8))
        ttk.Label(caja, text='Ojo: el sistema de mareas lo decide usted. Si el encabezado del archivo dice otra '
                             'cosa, el programa lo avisa en la bitácora pero usa lo que usted eligió. Para '
                             'EGM2008 descargado de ICGEM el encabezado indica «tide_free».',
                  style='Ayuda.TLabel', wraplength=1150).pack(anchor='w', pady=4)

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
        self.lbl_N.configure(text=explicar_N(m['sistema']))

    def _cambiar_sis_malla(self, _=None):
        s = self.tv_mallas.selection()
        if not s:
            messagebox.showinfo('Mallas', 'Primero seleccione una malla en la tabla.')
            return
        k = int(s[0])
        self.mallas[k]['sistema'] = self._sis(self.cb_sis_malla)
        self._refrescar_mallas()
        self.tv_mallas.selection_set(str(k))
        self.lbl_N.configure(text=explicar_N(self.mallas[k]['sistema']))

    # ---------------------------------------------------------- 4. Calcular
    def _tab_calcular(self):
        fr = ttk.Frame(self.nb)
        self.nb.add(fr, text='  4 · Calcular  ')
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
        filas = [('IDW: exponente p', 'idw_p', '0.5 – 6; mayor p = más peso a los nodos cercanos'),
                 ('IDW: número de nodos k', 'idw_k', '3 – 64'),
                 ('Kriging: número de vecinos n', 'krig_n', '4 – 64'),
                 ('Kriging: semiancho de ventana', 'krig_ventana', '3 – 12 nodos (para ajustar el variograma)'),
                 ('Depuración R6: umbral', 'umbral_mad', '× σ robusta; mayor = rechaza menos puntos')]
        for r, (t, k, ay) in enumerate(filas):
            ttk.Label(par, text=t).grid(row=r, column=0, sticky='w', padx=8, pady=3)
            ttk.Entry(par, textvariable=self.v[k], width=8).grid(row=r, column=1, sticky='w')
            ttk.Label(par, text=ay, style='Ayuda.TLabel').grid(row=r, column=2, sticky='w', padx=8)
        ttk.Label(par, text='Modelo paramétrico').grid(row=5, column=0, sticky='w', padx=8, pady=3)
        fm = ttk.Frame(par)
        fm.grid(row=5, column=1, columnspan=2, sticky='w')
        ttk.Radiobutton(fm, text='Ec. (18): δH + 2 inclinaciones (paper)', value='18',
                        variable=self.v['modelo']).pack(anchor='w')
        ttk.Radiobutton(fm, text='Ec. (23): superficie de 4 parámetros', value='23',
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
        self.txt_log = tk.Text(bit, wrap='none', font=('Consolas', 9), bg='#1b2631', fg='#e8eef4',
                               insertbackground=BLANCO, relief='flat', padx=8, pady=6)
        sy = ttk.Scrollbar(bit, command=self.txt_log.yview)
        sx = ttk.Scrollbar(bit, orient='horizontal', command=self.txt_log.xview)
        self.txt_log.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.pack(side='right', fill='y')
        sx.pack(side='bottom', fill='x')
        self.txt_log.pack(fill='both', expand=True)
        self.txt_log.tag_configure('mejor', foreground='#58d68d')
        self.txt_log.tag_configure('aviso', foreground='#f5b041')

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
                else:
                    self.btn_ejecutar.configure(state='normal')
                    self.btn_detener.configure(state='disabled')
                    mp.LOG_HOOK = None
                    if tipo == 'fin':
                        self.salida = msg[1]
                        self._mostrar_resultados()
                        self._preparar_kriging()
                        self.estado.set(f'Terminado. Resultados en {self.salida["carpeta"]}')
                        if self.salida['resultados']:
                            self.nb.select(4)
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
        self.nb.add(fr, text='  5 · Resultados  ')
        self._paso(fr, 'Resultados — modelo paramétrico y mejor método',
                   'Para cada malla: parámetros de la Ec. (18) con su desviación estándar, σ₀, R² y R² '
                   'ajustado de cada método de interpolación. En verde, el método con mayor R²_adj.')
        top = ttk.Frame(fr)
        top.pack(fill='x', padx=14)
        self.lbl_dep = ttk.Label(top, text='Aún no hay resultados: ejecute el cálculo en la pestaña 4.')
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
        t = tk.Text(f, height=12, font=('Consolas', 9), bg=BLANCO, relief='flat', padx=10, pady=8)
        t.pack(fill='both', expand=True)
        t.insert('end', 'Carpeta: ' + s['carpeta'] + '\n\n')
        expl = {'coordenadas_zero_tide.csv': 'puntos aceptados con h y H originales, corrección y valor zero tide',
                'puntos_rechazados.csv': 'puntos depurados y motivo (reglas R1-R7)',
                'resumen_modelos.csv': 'parámetros, σ0, R² y R²_adj de todos los métodos y mallas',
                'reporte.txt': 'bitácora completa'}
        for nombre in sorted(os.listdir(s['carpeta'])):
            e = expl.get(nombre, '')
            if nombre.startswith('interpolacion_'):
                e = 'N, l = h−H−N y residuos v de los 6 métodos en cada punto; σ y modelo del kriging'
            elif nombre.startswith('modelo_parametrico_'):
                e = 'tabla de la Ec. (18) para esa malla'
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
                  '  El modelo se ajustó por mínimos cuadrados con P = I (Ec. 19-20).',
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
            cv.create_text(x0 - 6, y + alto / 2, text=met, anchor='e', font=('Segoe UI', 8))
            cv.create_rectangle(x0, y + 3, x0 + largo, y + alto - 3, fill=color, outline='')
            cv.create_text(x0 + largo + 4, y + alto / 2, text=f'{v:.4f}', anchor='w',
                           font=('Segoe UI', 8, 'bold' if met == r['mejor'] else 'normal'))

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
        self.nb.add(fr, text='  6 · Kriging: semivariograma  ')
        self._paso(fr, 'Kriging ordinario — semivariograma / covariograma en un punto',
                   'Para cada punto el variograma se ajusta con los nodos de la malla en una ventana '
                   'alrededor del punto (2·semiancho × 2·semiancho nodos). Puntos = valores experimentales '
                   '(tamaño ∝ número de pares); líneas = los 4 modelos ajustados; línea gruesa = modelo elegido.')
        ctrl = ttk.Frame(fr)
        ctrl.pack(fill='x', padx=14)
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
        self.cb_k_punto.bind('<<ComboboxSelected>>', lambda e: self._dibujar_variograma())
        self.cb_k_punto.bind('<Return>', lambda e: self._dibujar_variograma())

        cuerpo = ttk.Frame(fr)
        cuerpo.pack(fill='both', expand=True, padx=14, pady=8)
        self.cv_var = tk.Canvas(cuerpo, bg=BLANCO, highlightthickness=1, highlightbackground='#d5dbe3')
        self.cv_var.pack(side='left', fill='both', expand=True)
        self.cv_var.bind('<Configure>', lambda e: self._dibujar_variograma())
        self.txt_var = tk.Text(cuerpo, width=46, font=('Consolas', 9), bg=BLANCO, relief='flat',
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

    def _datos_variograma(self):
        if not self.salida or not self.salida['resultados']:
            return None, 'Ejecute el cálculo (pestaña 4) para ver el semivariograma del kriging.'
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
                cv.create_text(X(v), H - B + 14, text=fmt_num(v), font=('Segoe UI', 8), fill=GRIS)
        for v in numeros_bonitos(ymin, ymax):
            if ymin <= v <= ymax:
                cv.create_line(L, Y(v), W - R, Y(v), fill='#eef1f5')
                cv.create_text(L - 6, Y(v), text=fmt_num(v), anchor='e', font=('Segoe UI', 8), fill=GRIS)
        cv.create_rectangle(L, T, W - R, H - B, outline='#9aa5b1')
        if ymin < 0:
            cv.create_line(L, Y(0), W - R, Y(0), fill='#9aa5b1', dash=(2, 2))
        cv.create_text((L + W - R) / 2, H - 18, text='distancia h  [km]  (proyección Origen-Nacional)',
                       font=('Segoe UI', 9), fill='#1b2631')
        cv.create_text(18, (T + H - B) / 2, text='C(h)  [m²]' if cov else 'γ(h)  [m²]', angle=90,
                       font=('Segoe UI', 9), fill='#1b2631')
        titulo = 'Covariograma' if cov else 'Semivariograma'
        cv.create_text(L, 22, anchor='w', fill=AZUL, font=('Segoe UI', 12, 'bold'),
                       text=f'{titulo} local — {r["malla"].nombre}, punto {pid}')
        # meseta (semivariograma) o C(0) (covariograma)
        elegido = next(t for t in mod['todos'] if t['tipo'] == mod['tipo'])
        meseta = elegido['nugget'] + elegido['sill']
        if ymin <= meseta <= ymax:
            cv.create_line(L, Y(meseta), W - R, Y(meseta), fill='#b03a2e', dash=(6, 4))
            cv.create_text(W - R - 4, Y(meseta) - 8 if not cov else Y(meseta) + 8, anchor='e',
                           fill='#b03a2e', font=('Segoe UI', 8),
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
        cv.create_text(lx + 16, ly, anchor='w', font=('Segoe UI', 8), text='experimental (tamaño ∝ nº pares)')
        for k, (t, _) in enumerate(curvas):
            yy = ly + 20 * (k + 1)
            es = t['tipo'] == mod['tipo']
            estilo = {'width': 3} if es else {'width': 1.5, 'dash': (5, 3)}
            cv.create_line(lx, yy, lx + 24, yy, fill=COLORES_MODELO[t['tipo']], **estilo)
            cv.create_text(lx + 32, yy, anchor='w', font=('Segoe UI', 8, 'bold' if es else 'normal'),
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
        for f in ('QGeoidCOL2023.txt', 'coords_EGM2008.txt'):
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
