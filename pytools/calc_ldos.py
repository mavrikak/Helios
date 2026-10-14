#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_ldos.py — Automated Dipole Decay Rate, Lifetime, and LDOS Calculator.

This script parses job XMLs to find dipole configurations, automatically
generates the required .pos files, triggers the HELIOS post-processing 
(smartly skipping if data already exists), and calculates emission properties.
"""

import argparse
import re
import csv
import sys
import subprocess
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

EPS_0 = 8.8541878128e-12  # F/m

def _parse_complex_triplet_line(line: str):
    s = line.strip()
    if not s or s[0] in "#%/": return None
    m_all = re.findall(r"\(\s*([-+]?[\d\.Ee+-]+)\s*,\s*([-+]?[\d\.Ee+-]+)\s*\)", s)
    if len(m_all) >= 3:
        return [complex(float(a), float(b)) for (a, b) in m_all[:3]]
    return None

def get_wavelengths_from_jobs(jobs_dir: Path) -> dict:
    lam_map = {}
    for xml_file in jobs_dir.glob("job.*.xml"):
        m_idx = re.search(r"job\.(\d+)\.xml", xml_file.name)
        if not m_idx: continue
        sim_idx = int(m_idx.group(1))
        
        txt = xml_file.read_text(errors="ignore")
        m_wav = re.search(r"<\s*wavelength\s*>([^<]+)<\s*/\s*wavelength\s*>", txt, re.IGNORECASE)
        
        if m_wav:
            try: 
                lam_map[sim_idx] = float(m_wav.group(1).strip())
            except ValueError: 
                pass
                
    if not lam_map:
        raise RuntimeError(f"No wavelengths could be parsed from XML files in {jobs_dir}")
        
    return lam_map

def get_all_dipole_subjobs(jobs_dir: Path) -> list:
    job_file = jobs_dir / "job.1.xml"
    if not job_file.exists():
        raise FileNotFoundError(f"Could not find {job_file}")
    
    txt = job_file.read_text(errors="ignore")
    subjobs = []
    
    block_m = re.search(r"<\s*job\b[\s\S]*?</\s*job\s*>", txt, re.IGNORECASE)
    if block_m and re.search(r"<\s*dipole\b", block_m.group(0), re.IGNORECASE):
        subjobs.append(1)
        
    for sm in re.finditer(r"<\s*subjob\b[\s\S]*?</\s*subjob\s*>", txt, re.IGNORECASE):
        sblk = sm.group(0)
        if re.search(r"<\s*dipole\b", sblk, re.IGNORECASE):
            lm = re.search(r"<\s*label\s*>\s*([^<]+)\s*</\s*label\s*>", sblk, re.IGNORECASE)
            if lm:
                m = re.search(r"inc_(\d+)", lm.group(1))
                if m: subjobs.append(int(m.group(1)))
                
    return sorted(list(set(subjobs)))

def get_dipoles_from_job(jobs_dir: Path, subjob_idx: int) -> list:
    job_file = jobs_dir / "job.1.xml"
    txt = job_file.read_text(errors="ignore")
    
    if subjob_idx == 1:
        block_m = re.search(r"<\s*job\b[\s\S]*?</\s*job\s*>", txt, re.IGNORECASE)
        block = block_m.group(0) if block_m else txt
    else:
        block = None
        for sm in re.finditer(r"<\s*subjob\b[\s\S]*?</\s*subjob\s*>", txt, re.IGNORECASE):
            sblk = sm.group(0)
            lm = re.search(r"<\s*label\s*>\s*([^<]+)\s*</\s*label\s*>", sblk, re.IGNORECASE)
            if lm and f"inc_{subjob_idx}" in lm.group(1):
                block = sblk
                break
                
    if not block:
        return []

    dipoles = []
    for dm in re.finditer(r"<\s*dipole\b[\s\S]*?</\s*dipole\s*>", block, re.IGNORECASE):
        dblk = dm.group(0)
        loc_m = re.search(r"<\s*location\s*>([^<]+)</\s*location\s*>", dblk, re.IGNORECASE)
        pol_m = re.search(r"<\s*polarization\s*>([^<]+)</\s*polarization\s*>", dblk, re.IGNORECASE)
        if loc_m and pol_m:
            loc = np.array([float(x) for x in loc_m.group(1).split()])
            pol = np.array([float(x) for x in pol_m.group(1).split()], dtype=np.complex128)
            dipoles.append({'loc': loc, 'pol': pol})
            
    return dipoles

def main():
    parser = argparse.ArgumentParser(description="Automated Dipole Decay Rate Pipeline.")
    parser.add_argument("sim", help="Simulation name (subfolder of sim_res)")
    parser.add_argument("--root", default="sim_res", help="Results root (default: sim_res)")
    parser.add_argument("--subjobs", type=str, default="all", help="Comma-separated list of subjobs (e.g., '2,3,4') or 'all' (default: all)")
    parser.add_argument("--n-bg", type=float, default=1.0, help="Background refractive index at dipole location (default: 1.0)")
    parser.add_argument("--lam-scale", type=float, default=1.0, help="Multiplier to convert wavelength to meters (e.g., 1e-9 if XML is in nm)")
    parser.add_argument("--points", type=str, default=None, help="Optional name for the generated points file")
    
    # Automation & Parallel arguments
    parser.add_argument("--parallel", action="store_true", help="Use dist_run_helios.py instead of run_sie.py")
    parser.add_argument("--mode", type=int, default=2, help="Simulation mode for dist_run_helios.py (default: 2)")
    parser.add_argument("--force-post", action="store_true", help="Force post-processing even if field files already exist")
    parser.add_argument("--skip-post", action="store_true", help="Strictly skip post-processing")
    
    # Output and Plotting arguments
    parser.add_argument("-o", "--output", type=str, default=None, help="Save results to CSV file")
    parser.add_argument("--plot", nargs="+", choices=["gamma", "tau", "ldos"], default=["gamma"], help="Quantities to plot.")
    parser.add_argument("--save", type=str, default=None, help="Save figure to this path")
    parser.add_argument("--show", action="store_true", help="Show interactive plot window")

    # Capture any unknown arguments to pass directly to the runner (like -th 8, -C 4, etc.)
    args, unknown_args = parser.parse_known_args()
    
    # Paths setup
    repo_root = Path(__file__).resolve().parent.parent
    root_dir = (repo_root / args.root).resolve()
    sim_dir = root_dir / args.sim
    jobs_dir = sim_dir / "jobs"
    points_dir = sim_dir / "points"
    
    run_sie_path = repo_root / "run_sie.py"
    dist_run_path = repo_root / "dist_run_helios.py"

    if not jobs_dir.exists():
        raise FileNotFoundError(f"Jobs directory not found: {jobs_dir}")

    lam_map = get_wavelengths_from_jobs(jobs_dir)

    if args.subjobs.lower() == "all":
        target_subjobs = get_all_dipole_subjobs(jobs_dir)
        print(f"[info] Auto-detected dipole subjobs: {target_subjobs}")
    else:
        target_subjobs = [int(x.strip()) for x in args.subjobs.split(",") if x.strip().isdigit()]

    global_metrics = []

    # Processing Loop
    for sj in target_subjobs:
        dipoles = get_dipoles_from_job(jobs_dir, sj)
        if not dipoles:
            print(f"[warn] No dipoles found in subjob {sj}. Skipping.")
            continue
            
        num_dipoles = len(dipoles)
        print(f"\n[info] --- Processing Subjob {sj} ({num_dipoles} dipoles) ---")
        
        # 1. Create .pos file automatically
        if args.points:
            pos_name = args.points
        else:
            pos_name = f"auto_dipoles_subjob_{sj}.pos"
            
        pos_path = points_dir / pos_name
        pos_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(pos_path, "w") as f:
            for dip in dipoles:
                f.write(f"{dip['loc'][0]:.6f} {dip['loc'][1]:.6f} {dip['loc'][2]:.6f}\n")
        print(f"[info] Points file generated at: {pos_path.relative_to(repo_root)}")

        # 2. Smart Post-Processing Logic
        fields_dir = sim_dir / "out" / "fields" / pos_path.stem
        existing_esc_files = list(fields_dir.glob(f"sim*_inc_{sj}.esc"))
        
        do_post = False
        if args.skip_post:
            print("[info] Skipping post-processing (--skip-post active).")
        elif args.force_post:
            print("[info] Forcing post-processing (--force-post active).")
            do_post = True
        elif not existing_esc_files:
            print("[info] No existing .esc files found. Initiating solver.")
            do_post = True
        else:
            print(f"[info] Found {len(existing_esc_files)} existing .esc files. Bypassing solver (use --force-post to override).")

        if do_post:
            if args.parallel:
                cmd = [
                    sys.executable, str(dist_run_path), 
                    "--mode", str(args.mode), 
                    "-s", args.sim, 
                    "-N", "n03,n04,n01", 
                    "-r", str(repo_root), 
                    "post", "-p", str(pos_path), 
                    "-a", "-th", "0"
                ] + unknown_args
            else:
                cmd = [
                    sys.executable, str(run_sie_path), 
                    "post", args.sim, 
                    "--mode", str(args.mode),
                    "-p", str(pos_path),
                    "-a", "-th", "0"
                ] + unknown_args
                
            print(f"[info] Running: {' '.join(cmd)}")
            subprocess.run(cmd, cwd=str(repo_root), check=True)

        # 3. Calculate metrics
        # Fetch the list of files again in case they were just generated
        esc_files = list(fields_dir.glob(f"sim*_inc_{sj}.esc"))
        
        if not esc_files:
            print(f"[warn] No .esc files available for subjob {sj}. Skipping calculations.")
            continue

        for f in esc_files:
            m = re.search(r"sim(\d+)_inc_", f.name)
            if not m: continue
            sim_idx = int(m.group(1))
            
            if sim_idx not in lam_map: continue
            # Convert wavelength to strict SI meters for EPS_0 math
            lam_val = lam_map[sim_idx] * args.lam_scale
            k_val = 2 * np.pi / lam_val
            
            E_scat_array = []
            with open(f, "r") as file:
                for line in file:
                    row = _parse_complex_triplet_line(line)
                    if row: E_scat_array.append(row)
                    
            if len(E_scat_array) < num_dipoles:
                continue
                
            E_scat_array = np.array(E_scat_array, dtype=np.complex128)

            total_im_dot = 0.0
            total_p_mag_sq = 0.0
            
            for i, dip in enumerate(dipoles):
                p_vec = dip['pol']
                E_scat = E_scat_array[i]
                
                p_mag_sq = np.sum(np.abs(p_vec)**2)
                im_dot = np.imag(np.vdot(p_vec, E_scat))
                
                total_im_dot += im_dot
                total_p_mag_sq += p_mag_sq
                
                prefactor = (6 * np.pi * EPS_0) / (args.n_bg * (k_val**3) * p_mag_sq)
                gamma_norm_i = 1.0 + (prefactor * im_dot)
                tau_norm_i = 1.0 / gamma_norm_i if gamma_norm_i != 0 else np.inf
                
                target_name = f"Subjob_{sj}_Dipole_{i+1}"
                global_metrics.append({
                    'sim_idx': sim_idx, 'lam': lam_val / args.lam_scale, 
                    'target': target_name, 'loc': dip['loc'], 
                    'im_dot': im_dot, 'gamma': gamma_norm_i, 'tau': tau_norm_i
                })

            if num_dipoles > 1:
                prefactor_tot = (6 * np.pi * EPS_0) / (args.n_bg * (k_val**3) * total_p_mag_sq)
                # WARNING: Assumes non-interacting dipoles in vacuum. Missing free-space mutual coupling.
                gamma_norm_tot = 1.0 + (prefactor_tot * total_im_dot)
                tau_norm_tot = 1.0 / gamma_norm_tot if gamma_norm_tot != 0 else np.inf
                
                global_metrics.append({
                    'sim_idx': sim_idx, 'lam': lam_val / args.lam_scale, 
                    'target': f"Subjob_{sj}_Coherent_Total", 'loc': np.array([np.nan, np.nan, np.nan]), 
                    'im_dot': total_im_dot, 'gamma': gamma_norm_tot, 'tau': tau_norm_tot
                })

    # --- GLOBAL MACROSCOPIC LDOS CALCULATION ---
    results = []
    from collections import defaultdict
    
    # Group globally by (Simulation Index, Rounded Location) across ALL subjobs
    loc_groups = defaultdict(list)
    for d in global_metrics:
        loc_key = (round(d['loc'][0], 4), round(d['loc'][1], 4), round(d['loc'][2], 4)) if not np.isnan(d['loc'][0]) else ("Total", d['target'])
        loc_groups[(d['sim_idx'], loc_key)].append(d)
        
    for key, group in loc_groups.items():
        # The macroscopic LDOS is the global average of decay rates at this location
        ldos_avg = np.mean([d['gamma'] for d in group])
        for d in group:
            results.append((d['sim_idx'], d['lam'], d['target'], d['im_dot'], d['gamma'], d['tau'], ldos_avg))

    if not results:
        print("\n[info] No data calculated. Exiting.")
        sys.exit(0)

    results.sort(key=lambda x: (x[1], x[2]))

    # Save CSV
    if args.output:
        out_path = Path(args.output)
        with open(out_path, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Sim_Index", "Wavelength", "Target", "Raw_Im_p_dot_E", "Gamma_norm", "Tau_norm", "LDOS_norm"])
            for row in results: writer.writerow(row)
        print(f"\n[saved] Metrics calculated and written to {out_path}")

    # Plotting
    plot_data = {}
    for row in results:
        target = row[2]
        if target not in plot_data:
            plot_data[target] = {'lam': [], 'gamma': [], 'tau': [], 'ldos': []}
        plot_data[target]['lam'].append(row[1])
        plot_data[target]['gamma'].append(row[4])
        plot_data[target]['tau'].append(row[5])
        plot_data[target]['ldos'].append(row[6])

    if plot_data and (args.show or args.save):
        plot_configs = []
        if "gamma" in args.plot:
            plot_configs.append(("gamma", r"Normalized Decay Rate ($\Gamma / \Gamma_0$)"))
        if "tau" in args.plot:
            plot_configs.append(("tau", r"Normalized Lifetime ($\tau / \tau_0$)"))
        if "ldos" in args.plot:
            plot_configs.append(("ldos", r"Normalized Proj. LDOS ($\rho / \rho_0$)"))

        num_plots = len(plot_configs)
        fig, axes = plt.subplots(num_plots, 1, figsize=(8.5, 3.8 * num_plots), sharex=True)
        if num_plots == 1: axes = [axes]

        for ax, (ptype, ylabel) in zip(axes, plot_configs):
            seen_ldos_signatures = set()
            for target, data in plot_data.items():
                sorted_pairs = sorted(zip(data['lam'], data[ptype]))
                x_vals, y_vals = zip(*sorted_pairs)
                if ptype == "ldos":
                    # Collapse duplicate LDOS curves into a single line
                    sig = tuple(np.round(y_vals, 4))
                    if sig in seen_ldos_signatures:
                        continue
                    seen_ldos_signatures.add(sig)
                    label_name = "Macroscopic LDOS (Avg)" if len(plot_data) > 1 else target
                    ax.plot(x_vals, y_vals, marker="o", label=label_name)
                else:
                    ax.plot(x_vals, y_vals, marker="o", label=target)
            
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.3)
            ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

        axes[-1].set_xlabel("Wavelength")
        axes[0].set_title("Automated Dipole Emission Properties")
        fig.tight_layout()
        
        if args.save:
            outp = Path(args.save)
            fig.savefig(outp, dpi=300, bbox_inches="tight")
            print(f"[saved] Figure written to {args.save}")
            
        if args.show:
            plt.show()
            
        plt.close(fig)

if __name__ == "__main__":
    main()