import os
import glob
import re

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: matplotlib is not installed. Please 'module load anaconda' first.")
    exit(1)

def extract_current_dir_data():
    """Finds all chunks in the current working directory and stitches them."""
    global_steps = []
    global_temps = []
    cumulative_step = 0
    
    # Dynamically find all OSZICAR files in chunk_* folders right here
    search_pattern = "chunk_*/OSZICAR"
    oszicar_files = glob.glob(search_pattern)
    
    # Sort them numerically (chunk_1, chunk_2, ... chunk_10)
    oszicar_files.sort(key=lambda x: int(re.search(r'chunk_(\d+)', x).group(1)))
    
    if not oszicar_files:
        print("[!] No chunk_*/OSZICAR files found in the current directory.")
        return global_steps, global_temps
        
    print(f"Scanning {len(oszicar_files)} chunks in {os.path.basename(os.getcwd())}...")
    
    for filepath in oszicar_files:
        chunk_step_count = 0
        if not os.path.getsize(filepath) > 0:
            continue
            
        with open(filepath, 'r') as f:
            for line in f:
                if "T=" in line and "E=" in line:
                    # Explicitly find the number immediately following 'T='
                    match = re.search(r'T=\s*([-\d\.]+)', line)
                    if match:
                        try:
                            temp = float(match.group(1))
                            # Ignore physically impossible temperatures from MPI text scrambling
                            if 0 < temp < 2000:
                                cumulative_step += 1
                                global_steps.append(cumulative_step)
                                global_temps.append(temp)
                                chunk_step_count += 1
                        except ValueError:
                            pass
        print(f"  -> {filepath}: Extracted {chunk_step_count} valid steps")
        
    return global_steps, global_temps

print("--- AIMD Local Trajectory Parser ---")
steps, temps = extract_current_dir_data()

if steps:
    # Build the plot
    plt.figure(figsize=(10, 6))
    
    # Convert to picoseconds (0.5 fs per step)
    ps = [s * 0.0005 for s in steps]
    
    plt.plot(ps, temps, label=f"Trajectory (Total: {ps[-1]:.2f} ps)", color="dodgerblue", linewidth=1.5)
    plt.axhline(y=298.15, color='black', linestyle='--', linewidth=2, label='Target (298.15 K)')
    
    # Formatting
    plt.xlabel("Simulated Time (Picoseconds)", fontsize=14, fontweight='bold')
    plt.ylabel("Temperature (K)", fontsize=14, fontweight='bold')
    plt.title(f"AIMD Continuous Heating Ramp\nDirectory: {os.path.basename(os.getcwd())}", fontsize=16, fontweight='bold')
    plt.legend(loc="lower right", fontsize=12)
    plt.grid(True, alpha=0.3)
    
    # Save the plot
    output_filename = "current_ramp.png"
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"------------------------------\nSuccess! Plot saved as {output_filename}")
