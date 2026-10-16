#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compara dos mallas de N nodo a nodo (A - B), por ejemplo la EGM2008 convertida
a zero tide por modelo_parametrico_zero_tide.py contra la descargada en zero
tide desde ICGEM.

Uso:
  python comparar_mallas.py resultados/coords_EGM2008_zero_tide.txt EGM2008_zero_tide_ICGEM.gdf
  python comparar_mallas.py A.txt B.txt --csv diferencias.csv
"""
import argparse
import csv
import math
import statistics

import modelo_parametrico_zero_tide as mp


def main():
    ap = argparse.ArgumentParser(description='Diferencias nodo a nodo entre dos mallas (A - B).')
    ap.add_argument('malla_a')
    ap.add_argument('malla_b')
    ap.add_argument('--csv', help='archivo CSV con lat, lon, N_A, N_B, A-B de cada nodo común')
    args = ap.parse_args()

    print('Leyendo mallas:')
    a = mp.cargar_malla(args.malla_a)
    b = mp.cargar_malla(args.malla_b)
    nodos_b = {(round(la, 6), round(lo, 6)): b.Z[i][j]
               for i, la in enumerate(b.lat) for j, lo in enumerate(b.lon) if b.Z[i][j] is not None}
    filas, dif = [], []
    for i, la in enumerate(a.lat):
        for j, lo in enumerate(a.lon):
            za = a.Z[i][j]
            zb = nodos_b.get((round(la, 6), round(lo, 6)))
            if za is None or zb is None:
                continue
            d = za - zb
            dif.append(d)
            filas.append((la, lo, za, zb, d))
    if not dif:
        raise SystemExit('Las mallas no tienen nodos en común (revise extensión y paso).')

    rms = math.sqrt(sum(d * d for d in dif) / len(dif))
    k = max(range(len(dif)), key=lambda n: abs(dif[n]))
    print(f'\nDiferencia A - B en {len(dif)} nodos comunes '
          f'(A: {len(a.lat) * len(a.lon)}, B: {len(b.lat) * len(b.lon)} nodos):')
    print(f'  media  {statistics.fmean(dif) * 1000:9.3f} mm')
    print(f'  desv.  {statistics.pstdev(dif) * 1000:9.3f} mm')
    print(f'  RMS    {rms * 1000:9.3f} mm')
    print(f'  mín    {min(dif) * 1000:9.3f} mm')
    print(f'  máx    {max(dif) * 1000:9.3f} mm')
    print(f'  |máx| en lat {filas[k][0]:.4f}, lon {filas[k][1]:.4f}')

    # Diferencia media por banda de latitud (la marea depende solo de φ)
    bandas = {}
    for la, _, _, _, d in filas:
        bandas.setdefault(math.floor(la), []).append(d)
    print('\n  Media de A - B por banda de latitud de 1°:')
    for g in sorted(bandas):
        print(f'    [{g:>4}°, {g + 1:>4}°)  {statistics.fmean(bandas[g]) * 1000:9.3f} mm  '
              f'({len(bandas[g])} nodos)')

    if args.csv:
        with open(args.csv, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['lat', 'lon', 'N_A', 'N_B', 'A_menos_B'])
            for la, lo, za, zb, d in filas:
                w.writerow([f'{la:.6f}', f'{lo:.6f}', f'{za:.6f}', f'{zb:.6f}', f'{d:.6f}'])
        print(f'\nDiferencias por nodo guardadas en {args.csv}')


if __name__ == '__main__':
    main()
