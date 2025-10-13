#!/usr/bin/env python3
"""
gen_vienna_hairpins.py v2.6

Genera M secuencias P(14)-Q(4)-R(14)-S(21), valida con ViennaRNA (RNAfold),
rechaza secuencias que formen más de 1 horquilla o donde Q empareje con P/R,
genera imágenes (RNAplot + gs), guarda resultados en formato Vienna.

Compatible con ViennaRNA 2.7.0

Requisitos:
 - RNAfold y RNAplot en PATH (ViennaRNA 2.7+)
 - Ghostscript 'gs' (para PNG)
"""

import random
import subprocess
import shutil
import os
import sys
import glob
from typing import Set, Tuple, List, Dict

# Parámetros del esquema (fijos según tu especificación)
P = 14
Q = 4
R = 14
S = 21
MIN_LEN = P + Q + R + S  # 53

# Random seed para reproducibilidad (opcional)
random.seed()

# Complemento ADN
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def complement_base(b: str) -> str:
    return COMP.get(b.upper(), 'N')

def complement_reverse(seq: str) -> str:
    return ''.join(complement_base(c) for c in reversed(seq))

def random_seq(length: int) -> str:
    return ''.join(random.choice(['A','C','G','T']) for _ in range(length))

# --- Dot-bracket -> pares (i,j) ---
def pairs_from_dotbracket(db: str) -> Set[Tuple[int,int]]:
    stack = []
    pairs = set()
    for i, ch in enumerate(db):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if not stack:
                continue
            j = stack.pop()
            a, b = (j, i) if j < i else (i, j)
            pairs.add((a,b))
    return pairs

# --- Pares diseñados P<->R (índices 0-based) ---
def designed_pairs_indices(P_len: int, Q_len: int, R_len: int) -> Set[Tuple[int,int]]:
    """
    Para i in [0..P_len-1] -> j = P_len + Q_len + (R_len-1 - i)
    """
    pairs = set()
    for i in range(min(P_len, R_len)):
        j = P_len + Q_len + (R_len - 1 - i)
        a, b = (i, j) if i < j else (j, i)
        pairs.add((a,b))
    return pairs

DESIGNED_PAIRS = designed_pairs_indices(P, Q, R)

# --- Ejecutar RNAfold --noPS sobre una secuencia y parsear salida ---
def evaluate_with_rnafold(seq: str) -> Tuple[str, float]:
    if not shutil.which('RNAfold'):
        raise RuntimeError("RNAfold no encontrado en PATH.")
    try:
        proc = subprocess.run(
            ['RNAfold', '--noPS'], 
            input=seq.encode(), 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            check=True
        )
        out = proc.stdout.decode().strip().splitlines()
        if len(out) < 2:
            raise RuntimeError("Salida inesperada de RNAfold")
        
        struct_line = out[1].strip()
        idx = struct_line.rfind('(')
        if idx == -1:
            parts = struct_line.split()
            struct = parts[0]
            energy = float(parts[-1].strip('()'))
        else:
            struct = struct_line[:idx].strip()
            energy = float(struct_line[idx:].strip('() ').replace(',', '.'))
        return struct, energy
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error RNAfold: {e.stderr.decode()}")

# --- Generar imagen usando RNAplot (v2.7 compatible) y convertir a PNG ---
def plot_structure(seq: str, struct: str, seq_idx: int) -> Tuple[str, str]:
    """
    Genera archivos seq_{idx}_ss.ps y seq_{idx}.png en carpeta plots/
    Compatible con RNAplot 2.7.0 (sin opciones -o)
    Retorna (ps_path, png_path)
    """
    if not shutil.which('RNAplot'):
        raise RuntimeError("RNAplot no encontrado en PATH.")
    
    # Verificar que longitudes coincidan
    if len(seq) != len(struct):
        raise RuntimeError(f"Longitudes no coinciden: seq={len(seq)}, struct={len(struct)}")
    
    # Nombres de archivos finales
    base = f"seq_{seq_idx}"
    ps_path = os.path.join('plots', f"{base}_ss.ps")
    png_path = os.path.join('plots', f"{base}.png")
    
    input_text = f"{seq}\n{struct}\n"
    
    # Guardar directorio original
    original_dir = os.getcwd()
    
    try:
        # Cambiar a directorio plots para que RNAplot genere archivos ahí
        os.chdir('plots')
        
        # Limpiar archivos PS/EPS antiguos en el directorio
        for old_file in glob.glob('*.ps') + glob.glob('*.eps'):
            try:
                os.remove(old_file)
            except:
                pass
        
        # RNAplot 2.7 lee stdin y genera rna_ss.ps
        proc = subprocess.run(
            ["RNAplot"],
            input=input_text,
            text=True,
            capture_output=True
        )
        
        # Verificar si hubo error
        if proc.returncode != 0:
            raise RuntimeError(f"RNAplot falló: {proc.stderr}")
        
        # Esperar un momento para que el archivo se escriba completamente
        import time
        time.sleep(0.1)
        
        # Buscar archivo .ps o .eps generado (intentar varios nombres comunes)
        possible_names = ['rna_ss.ps', 'rna.ps', 'rna.eps', 'ss.ps']
        generated_ps = None
        
        for name in possible_names:
            if os.path.exists(name):
                generated_ps = name
                break
        
        # Si no encontramos ninguno, buscar cualquier .ps/.eps
        if not generated_ps:
            ps_files = glob.glob('*.ps') + glob.glob('*.eps')
            if ps_files:
                generated_ps = ps_files[0]
        
        if not generated_ps:
            # Debug: listar todos los archivos
            all_files = os.listdir('.')
            raise RuntimeError(f"RNAplot no generó archivo PS/EPS. Archivos en plots/: {all_files}. stderr: {proc.stderr}")
        
        # Archivo encontrado
        final_ps_name = f"{base}_ss.ps"
        
        # Renombrar si es necesario
        if generated_ps != final_ps_name:
            if os.path.exists(final_ps_name):
                os.remove(final_ps_name)
            os.rename(generated_ps, final_ps_name)
        
        ps_path = os.path.join(os.getcwd(), final_ps_name)
        
        # Convertir a PNG con Ghostscript
        png_created = ''
        if shutil.which('gs'):
            try:
                png_name = f"{base}.png"
                # Ghostscript necesita -dEPSCrop para archivos EPS
                result = subprocess.run([
                    "gs", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET",
                    "-dEPSCrop",  # Importante para archivos EPS
                    "-sDEVICE=pngalpha", "-r300",  # Aumentar resolución a 300 DPI
                    f"-sOutputFile={png_name}", final_ps_name
                ], capture_output=True)
                
                if result.returncode != 0:
                    raise RuntimeError(f"Ghostscript error: {result.stderr.decode()}")
                
                if os.path.exists(png_name):
                    png_created = os.path.join(os.getcwd(), png_name)
                else:
                    raise RuntimeError(f"PNG no se generó: {png_name}")
                    
            except Exception as e:
                print(f"  ⚠️  Ghostscript falló para seq_{seq_idx}: {e}")
        else:
            print(f"  ⚠️  'gs' no encontrado; solo .ps disponible para seq_{seq_idx}")
        
        return ps_path, png_created
        
    except Exception as e:
        raise RuntimeError(f"Error en plot_structure: {e}")
    finally:
        # Volver al directorio original
        os.chdir(original_dir)

# --- Validación de estructura ---
def analyze_structure_validity(struct: str, seq_len: int, P_len: int, Q_len: int, R_len: int) -> Tuple[bool, str]:
    """
    Retorna (valid, message). valid=True si:
      - Los pares coinciden exactamente con DESIGNED_PAIRS
      - Q no forma pares con P o R
    """
    pred_pairs = pairs_from_dotbracket(struct)
    extra_pairs = pred_pairs - DESIGNED_PAIRS
    missing_pairs = DESIGNED_PAIRS - pred_pairs

    q_start = P_len
    q_end = P_len + Q_len - 1
    q_pairs = [p for p in pred_pairs if (q_start <= p[0] <= q_end) or (q_start <= p[1] <= q_end)]

    if extra_pairs:
        return False, f"Pares extra detectados ({len(extra_pairs)})."
    if missing_pairs:
        return False, f"Pares faltantes del stem ({len(missing_pairs)})."
    if q_pairs:
        return False, f"Loop Q presenta emparejamientos no permitidos."
    return True, "Válida"

# --- Generador principal ---
def generate_M_sequences(M: int, max_attempts=10000) -> List[Dict]:
    """
    Genera hasta M secuencias válidas.
    Retorna lista de dicts: id, seq, struct, energy, ps, png
    """
    if M <= 0:
        raise ValueError("M debe ser > 0")

    os.makedirs("plots", exist_ok=True)
    accepted = []
    seen = set()
    attempts = 0

    while len(accepted) < M and attempts < max_attempts:
        attempts += 1
        
        # Generar secuencia P-Q-R-S
        p_seq = random_seq(P)
        q_seq = random_seq(Q) if Q > 0 else ''
        r_seq = complement_reverse(p_seq)
        s_seq = random_seq(S) if S > 0 else ''
        full_seq = p_seq + q_seq + r_seq + s_seq
        
        if full_seq in seen:
            continue
        seen.add(full_seq)

        # Evaluar con RNAfold
        try:
            struct, energy = evaluate_with_rnafold(full_seq)
        except Exception as e:
            print(f"❌ Error RNAfold: {e}")
            continue

        # Validar estructura
        valid, msg = analyze_structure_validity(struct, len(full_seq), P, Q, R)
        if not valid:
            print(f"[RECHAZADA] intento {attempts}: {msg} | {p_seq}... | E={energy:.2f}")
            continue

        # Generar imágenes
        seq_id = len(accepted) + 1
        ps_path, png_path = '', ''
        try:
            ps_path, png_path = plot_structure(full_seq, struct, seq_id)
            print(f"[✓ {len(accepted)+1}/{M}] ID={seq_id} E={energy:.2f} {full_seq[:12]}... PNG={'✓' if png_path else '✗'}")
        except Exception as e:
            print(f"[✓ {len(accepted)+1}/{M}] ID={seq_id} E={energy:.2f} {full_seq[:12]}... PNG=✗ (Error: {e})")

        rec = {
            'id': seq_id,
            'seq': full_seq,
            'struct': struct,
            'energy': energy,
            'ps': ps_path,
            'png': png_path
        }
        accepted.append(rec)

    if len(accepted) < M:
        print(f"\n⚠️  Advertencia: {len(accepted)}/{M} secuencias generadas en {attempts} intentos.")
    else:
        print(f"\n✅ {len(accepted)} secuencias válidas en {attempts} intentos.")
    return accepted

# --- Guardar archivos ---
def save_vienna(filename: str, records: List[Dict]):
    with open(filename, 'w') as f:
        for r in records:
            f.write(r['seq'] + '\n')
            f.write(f"{r['struct']} ({r['energy']:.2f})\n")
    print(f"✅ Vienna: {filename}")

def save_csv_summary(filename: str, records: List[Dict]):
    with open(filename, 'w', newline='') as f:
        f.write("id,seq,struct,energy,ps,png\n")
        for r in records:
            f.write(f"{r['id']},{r['seq']},{r['struct']},{r['energy']},{r['ps']},{r['png']}\n")
    print(f"✅ CSV: {filename}")

def select_most_stable(records: List[Dict], top_k: int = 5):
    return sorted(records, key=lambda x: x['energy'])[:top_k]

# --- CLI principal ---
def main():
    print("=" * 70)
    print("  Generador P(14)-Q(4)-R(14)-S(21) + Validación ViennaRNA 2.7")
    print("=" * 70)
    
    try:
        M = int(input("\nIngresa M (número de secuencias a generar): ").strip())
    except:
        print("❌ Entrada inválida.")
        return

    # Verificar herramientas
    missing = []
    if not shutil.which('RNAfold'): missing.append('RNAfold')
    if not shutil.which('RNAplot'): missing.append('RNAplot')
    if missing:
        print(f"❌ Herramientas faltantes: {', '.join(missing)}")
        return
    
    if not shutil.which('gs'):
        print("⚠️  Ghostscript (gs) no encontrado. Solo se generarán .ps\n")
    
    print(f"\n🔬 Generando {M} secuencias válidas (longitud={MIN_LEN})...\n")
    
    records = generate_M_sequences(M)
    if not records:
        print("❌ No se generaron secuencias válidas.")
        return

    # Guardar resultados
    save_vienna('resultados.vienna', records)
    save_csv_summary('resultados_summary.csv', records)

    # Ranking por estabilidad
    top = select_most_stable(records, top_k=min(5, len(records)))
    print("\n" + "=" * 70)
    print("🏆 TOP ESTRUCTURAS MÁS ESTABLES (energía libre de Gibbs)")
    print("=" * 70)
    for i, t in enumerate(top, 1):
        png_status = "✓" if t['png'] else "✗"
        print(f"{i}. ID={t['id']:3d} | E={t['energy']:6.2f} kcal/mol | {t['seq'][:15]}... | PNG={png_status}")
    
    png_count = sum(1 for r in records if r['png'])
    print(f"\n📊 Resumen: {len(records)} secuencias | {png_count} imágenes PNG generadas")
    print(f"📁 Archivos en: plots/ | resultados.vienna | resultados_summary.csv")
    print("=" * 70)

if __name__ == "__main__":
    main()
