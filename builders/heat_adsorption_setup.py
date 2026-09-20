import os
import sys

# The HEAT package used by the original workflow is external to this repository.
# Set HEAT_APPS to the directory that makes `import heat` available if needed.
heat_apps = os.environ.get("HEAT_APPS")
if heat_apps:
    sys.path.append(heat_apps)
from heat import heat_surface

from ase.io import read

##### HEAT ADSORPTION CODE ######
#### INPUTS ####
#Path to clean slab calculation
clean_slab_path = "clean_slab" # Path to folder containing clean slab calculation

# Path to adsorbate POSCAR
adsorbate_poscar_path = "ads.POSCAR" # Path to POSCAR containing only the adsorbate molecule 

# Adsorbed atom index
adsorbed_atom_index = 0 # Index of atom in adsorbate POSCAR that forms bonds with surface, count starts from 0

# Types of sites to find (list containing one or more of "ontop", "bridge", and "hollow")
sites = ["ontop","bridge","hollow"]

# Number of adsorbates (list containing numbers of adsorbates) (for e.g., specifying 1 and 2 will
# lead to creation of structures with both 1 and 2 adsorbates on the surface)
n_adsorbates = [1]

# z-height difference between topmost atom in top layer and bottom most atom in top layer
# This is 0.9 for regular facets but needs to be changed to 2.1 for steps
diff_height = 1.1

# Whether spin polarization is enabled (set to True if ISPIN=2 and also specify magmom below)
spin_pol=False

# Whether to submit calculations after submitting? (recommended: False)
# Only set to True if you know all structures generated are going to be correct
submit = False

#### VASP inputs ####
# Dictionary containing magmom guesses for INCAR (element: magmom guess) (Keep this empty if ISPIN is 1)
#magmom={"Pt":0.01,"Fe":4,"Ni":2,"Co":3,"Cu":0.01,"Ir":0.01,"O":1.0,"H":0.01}
magmom={}

# Dictionary containing inputs for submission file 
run= {\
"QUEUE" : "standby" ,\
"WALLTIME" : "04:00:00" ,\
"RESUBMISSIONS" : "5" ,\
"NODES": "1" ,\
"CORES":"64",\
"CLUSTER": os.environ.get("HEAT_CLUSTER", "cluster"),\
"EMAIL": os.environ.get("HEAT_EMAIL", "user@example.edu")}

## Dictionary containing inputs for INCAR (add or remove tags as necessary)
incar= {\
"ENCUT" : "400" ,\
"EDIFF" : "1e-5"  ,\
"EDIFFG" : "-0.02" ,\
"IBRION" : "2" ,\
"GGA" : "PE" ,\
"ISPIN" : "1"  ,\
#"MAGMOM" : "0" ,\
"ISMEAR" : "1" ,\
"SIGMA" : "0.2" ,\
"ISIF" : "2" ,\
#"POTIM" : "0.25" ,\
#"LORBIT" : "11",\
#"ISTART" : "0" ,\
#"ICHARG" : "2" ,\
"ALGO" : "Normal",\
#"NWRITE" : "2",\
"ISYM" : "0" ,\
"PREC": "Normal",\
"NSW" : "500" ,\
"LCHARG" : "F" ,\
"LWAVE" : "F" ,\
"LREAL" : "F" ,\
"NELMIN" : "5" ,\
}

# Dictionary containing inputs for KPOINTS (add or remove tags as necessary)
kpoints={\
"KPOINTS" : (4,4,1) ,\
"TYPE" : "g" ,\
"SHIFT" : (0,0,0) ,\
}


##################### END INPUTS ###############################
################# DO NOT TOUCH ANYTHING BELOW ##################
# Read adsorbate
adsorbate = read(adsorbate_poscar_path)


#Make calculations
hea=heat_surface.HEASurface(incar=incar,kpoints=kpoints,run=run,magnetic_moments=magmom)
hea.create_adsorption_calculations_general(clean_slab_path,adsorbate=adsorbate,adsorbed_atom=adsorbed_atom_index,n_adsorbates=n_adsorbates,spin_pol=spin_pol,sites=sites,height=diff_height,submit=submit)

#### Go to the clean slab path to find the adsorption calculations
