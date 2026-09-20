import matplotlib.pyplot as plt
import os

print("Loading diagnostic data...")
if not os.path.exists('aimd_diagnostics.dat'):
    print("Error: aimd_diagnostics.dat not found. Run the grep command first.")
    exit(1)

steps, temp, energy = [], [], []

# Safely parse the file line-by-line
with open('aimd_diagnostics.dat', 'r') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 3:
            try:
                s = float(parts[0])
                t = float(parts[1])
                e = float(parts[2])  # Column 3 is the Total Energy
                steps.append(s)
                temp.append(t)
                energy.append(e)
            except ValueError:
                continue

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
