import os
import glob
import re

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: matplotlib is not installed. Please 'module load anaconda' first.")
    exit(1)

def get_target_steps():
    """Extracts TARGET_STEPS from submit.sh"""
    try:
        with open("submit.sh", "r") as f:
            for line in f:
                if line.startswith("TARGET_STEPS="):
                    return int(line.split("=")[1].strip())
    except Exception:
        pass
    return 5000

def extract_all_indicators():
    global_steps, temps, free_energies, kinetic_energies = [], [], [], []
    cumulative_step = 0
    
    search_pattern = "chunk_*/OSZICAR"
    oszicar_files = glob.glob(search_pattern)
    oszicar_files.sort(key=lambda x: int(re.search(r'chunk_(\d+)', x).group(1)))
    
    print(f"Scanning {len(oszicar_files)} chunks with MPI-Duplicate Filtering...")
    
    for filepath in oszicar_files:
        if not os.path.getsize(filepath) > 0:
            continue
            
        seen_steps_in_chunk = set() # Tracks unique steps to ignore MPI spam
        
        with open(filepath, 'r') as f:
            for line in f:
                if "T=" in line and "F=" in line and "EK=" in line:
                    # Capture the ionic step number at the start of the line
                    m_step = re.search(r'^\s*(\d+)\s+T=', line)
                    m_T = re.search(r'T=\s*([-\d\.]+)', line)
                    m_F = re.search(r'F=\s*([-\d\.\+E]+)', line)
                    m_EK = re.search(r'EK=\s*([-\d\.\+E]+)', line)
                    
                    if m_step and m_T and m_F and m_EK:
                        ionic_step = int(m_step.group(1))
                        
                        # Only process if we haven't seen this step number in this chunk yet
                        if ionic_step not in seen_steps_in_chunk:
                            try:
                                t_val = float(m_T.group(1))
                                f_val = float(m_F.group(1))
                                ek_val = float(m_EK.group(1))
                                
                                if 0 < t_val < 2000 and f_val < 0:
                                    seen_steps_in_chunk.add(ionic_step)
                                    cumulative_step += 1
                                    
                                    global_steps.append(cumulative_step)
                                    temps.append(t_val)
                                    free_energies.append(f_val)
                                    kinetic_energies.append(ek_val)
                            except ValueError:
                                pass
                                
    return global_steps, temps, free_energies, kinetic_energies

# Run extraction
target_steps = get_target_steps()
steps, temps, free_energies, kinetic_energies = extract_all_indicators()

if steps:
    total_steps_done = steps[-1]
    total_ps = total_steps_done * 0.0005
    target_ps = target_steps * 0.0005
    percent_complete = (total_steps_done / target_steps) * 100
    
    print("-" * 40)
    print(f"Target:     {target_steps} steps ({target_ps:.2f} ps)")
    print(f"TRUE Steps: {total_steps_done} steps ({total_ps:.4f} ps)")
    print(f"Progress:   {percent_complete:.1f}%")
    print("-" * 40)

    ps = [s * 0.0005 for s in steps]
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    
    fig.suptitle(f"AIMD Diagnostics: {os.path.basename(os.getcwd())}\n"
                 f"True Progress: {total_steps_done} / {target_steps} Steps ({percent_complete:.1f}%) | Time: {total_ps:.2f} ps", 
                 fontsize=16, fontweight='bold')
    
    ax1.plot(ps, temps, color="dodgerblue", linewidth=1.5)
    ax1.axhline(y=298.15, color='black', linestyle='--', label='Target (298 K)')
    ax1.set_ylabel("Temperature (K)", fontweight='bold')
    ax1.legend(loc="lower right")
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(ps, kinetic_energies, color="orange", linewidth=1.5)
    ax2.set_ylabel("Kinetic Energy (eV)", fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    ax3.plot(ps, free_energies, color="crimson", linewidth=1.5)
    ax3.set_ylabel("Free Energy F (eV)", fontweight='bold')
    ax3.set_xlabel("Simulated Time (Picoseconds)", fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("diagnostics.png", dpi=300)
    print("Success! Plot saved as diagnostics.png")
