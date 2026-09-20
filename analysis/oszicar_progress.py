import matplotlib.pyplot as plt
import glob
import os

print("Scanning OSZICAR files to extract diagnostic data...")

# Find all chunk directories and sort them mathematically
chunks = sorted(glob.glob('chunk_*'), key=lambda x: int(x.split('_')[1]))

if not chunks:
    print("Error: No chunk directories found.")
    exit(1)

steps, temp, energy = [], [], []
step_counter = 1

for chunk in chunks:
    oszicar_path = os.path.join(chunk, 'OSZICAR')
    if not os.path.exists(oszicar_path):
        continue
        
    with open(oszicar_path, 'r') as f:
        for line in f:
            # Look for the temperature line indicator
            if " T= " in line:
                parts = line.split()
                # Ensure the line isn't malformed and matches the expected layout
                if len(parts) >= 5 and parts[1] == "T=":
                    try:
                        t = float(parts[2])  # Temperature
                        e = float(parts[4])  # Total Energy
                        
                        # Filter out wild unphysical spikes from crash frames
                        if 0 < t < 2000:
                            steps.append(step_counter)
                            temp.append(t)
                            energy.append(e)
                            step_counter += 1
                    except ValueError:
                        continue # Skip glitched VASP formatting strings

if not steps:
    print("Error: No valid temperature data found in the OSZICAR files.")
    exit(1)

print(f"Successfully extracted {len(steps)} data points. Generating plot...")

# Create a stacked graph with 2 subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# --- Top Panel: Temperature ---
ax1.plot(steps, temp, color='#007acc', linewidth=1.2, alpha=0.9)
ax1.axhline(y=298.15, color='black', linestyle='--', linewidth=2, label='Target (298.15 K)')
ax1.set_title('AIMD Equilibration Diagnostics', fontsize=14, fontweight='bold')
ax1.set_ylabel('Temperature (K)', fontsize=12)
ax1.legend(loc='upper right')
ax1.grid(True, linestyle=':', alpha=0.6)

# --- Bottom Panel: Total Energy ---
ax2.plot(steps, energy, color='#e67e22', linewidth=1.2, alpha=0.9)
ax2.set_xlabel('Cumulative MD Steps', fontsize=12)
ax2.set_ylabel('Total Energy (eV)', fontsize=12)
ax2.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('diagnostics_full.png', dpi=300)
print("Success! High-resolution double graph saved as 'diagnostics_full.png'.")
