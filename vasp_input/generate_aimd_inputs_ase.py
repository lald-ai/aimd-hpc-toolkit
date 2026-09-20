from ase.calculators.vasp import Vasp
from ase.io import read, write

def main():
    infile = 'restart.json'
    ibrion = 0  # 2=ionic relaxation; 5=vibrations; -1=Single Point; 0=AIMD
    ncore = 16  # number of compute cores that work on an individual orbital
    npar = 8  # determines the number of bands that are treated in parallel.
    kpar = 5  # determines the number of k-points that are to be treated in parallel.

    # Electronic Convergence Criteria
    ediff = 1E-04  # Convergence for electronic SC loop.
    nelm = 200  # maximum number of electronic SCF steps.
    nelmin = 4  # minimum number of electronic SCF steps.
    # algo = 'Normal'  # electronic minimization algorithm
    ialgo = 48  # algorithm to optimize the orbitals
    prec = 'Normal'  # precision accuracy.
    lreal = 'Auto'  # Recommendation: False for bulk; Auto for surface.
    ispin = 1  # 1=non-spin polarized; 2=spin-polarized
    # idipol = 3  # 1-3: dipole moment will be calculated parallel to a,b,c respectively. 4: full dipole moment in all directions.

    # Ionic Relaxation Convergence Criteria
    ediffg = -0.1  # convergence criteria. +ve: energy tolerance; -ve: force tolerance.
    nsw = 875  # maximum number of ionic steps. or number of steps in MD runs
    isif = 0  # 2=cell shape and volume fixed, while positions are free to relax. 0=MD (for IBRION=0)

    # Molecular Dynamics Settings
    mdalgo = 0  # molecular-dynamics-simulation protocol. 2: Nose-Hoover thermostat
    potim = 1  # time step in fs for molecular dynamics
    tebeg = 0  # starting temperature in K
    teend = 350  # final temperature in K
    smass = -1  # controls the velocities during an AIMD run. 
    
    # VASP Calculator Parameters
    encut = 400  # plane-wave cut-off in eV
    xc = 'BEEF-vdW'  # xc-functional to be used
    # aggac = 0.0  # multiplier to gradient correction in the GGA correlation functional
    luse_vdw = True  # nonlocal vdW-DF functional
    lasph = True  # non-spherical contributions related to the gradient of the density in the PAW spheres.
    isym = 0  # Symmetry. 0:No Symmetry
    kpts = (4, 4, 1)  # k-point grid sub-divisions, k-point grid density
    gamma = True  # Gamma-point centered k-point sampling (defaults to Monkhorst-Pack)
    ismear = 1  # 1=MP order 1 (metals), 2=MP order 2, -1=Fermi, 0=Gaussian
    sigma = 0.2  # width of smearing. Default value of 0.2 is good for metals.
    pbc = (True, True, True)
    icharg = 2  # 1: Read the charge density from CHGCAR file, 2: superposition of atomic charge densities.
    lwave = False  # whether to write wavefunction to WAVECAR
    lcharg = False  # whether to write charge densities to CHGCAR
    nwrite = 1  # tag determines how much will be written to the file OUTCAR
    lblueout = True  # whether to write output for the free-energy gradient calculation to the REPORT file

    # Set up VASP calculator
    calc = Vasp(encut=encut,
                ibrion=ibrion,
                ncore=ncore,
                npar=npar,
                kpar=kpar,
                ediff=ediff,
                nelm=nelm,
                nelmin=nelmin,
                # algo=algo,
                ialgo=ialgo,
                prec=prec,
                lreal=lreal,
                ispin=ispin,
                # idipol=idipol,
                ediffg=ediffg,
                nsw=nsw,
                isif=isif,
                xc=xc,
                # aggac=aggac,
                luse_vdw=luse_vdw,
                lasph=lasph,
                isym=isym,
                kpts=kpts,
                gamma=gamma,
                ismear=ismear,
                sigma=sigma,
                mdalgo=mdalgo,
                potim=potim,
                tebeg=tebeg,
                teend=teend,
                smass=smass,
                icharg=icharg,
                lwave=lwave,
                lcharg=lcharg,
                nwrite=nwrite,
                lblueout=lblueout,
                )

    # Read in system coordinates
    atoms = read(infile)
    atoms.set_pbc(pbc)

    # Set calculator and perform calculation
    atoms.set_calculator(calc)
    calc.write_input(atoms)
    return None

if __name__ == '__main__':
    main()
