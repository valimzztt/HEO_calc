from ase.calculators.espresso import Espresso
from ase.db import connect
from clease.tools import update_db
from clease.settings import Concentration
import matplotlib.pyplot as plt
import json
conc = Concentration(basis_elements=[['Ti', 'Zr'], ['O']])
conc.set_conc_ranges(ranges=[[(0.5,0.5),(0.5,0.5)], [(1,1)]])

#define crystal structure
from clease.settings import CECrystal
settings = CECrystal(concentration=conc,
    spacegroup=60,
    basis=[(0.00000, 0.67652, 0.25000), (0.22868, 0.61936, -0.08125)],
    cell=[4.8127, 5.8570, 5.1086, 90, 90, 90],
    supercell_factor=2,
    db_name="clease.db",
    basis_func_type='binary_linear',
    #max_cluster_size=2,
    max_cluster_dia=[4])
pseudo_dir = './pseudos/'



pseudopotentials = {'Ti': pseudo_dir + 'Ti.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                    'Zr': pseudo_dir + 'Zr.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                    'O': pseudo_dir + 'O.pbesol-n-kjpaw_psl.1.0.0.UPF'}

input_data = {
        'control': {
           'calculation': 'relax',
           'restart_mode': 'from_scratch',
           'prefix': 'tutorial',     
            },
        'system': {
           'ecutwfc':  40.424222147,
           'occcupations': 'smearing',
           'degauss': 0.0146997171,
            }
}


calc = Espresso(pseudopotentials = pseudopotentials,
                input_data = input_data,
                kspacing = 0.4)

# Misc text
#    lreal=False,
#    lscalapack=False,
#    kspacing=0.4,
#    ediffg=-0.02,
#    encut=550,
#    ibrion=2,
#    isif=4,
#    nsw=150,
#    lasph=True,
#    ncore=4,
#    istart=0,
#    isym=0,
#    lwave=False,
#    ldau_luj={'Ti': {'L': 2, 'U': 5,'J': 0}, 'Zr': {'L':  2, 'U': 3, 'J': 0}})

#Setting up the database
db_name = "clease.db"
db = connect(db_name)

# Run calculations for all structures that are not converged.
for row in db.select(converged=False):
  atoms = row.toatoms()
  atoms.calc = calc
  atoms.get_potential_energy()
  update_db(uid_initial=row.id, final_struct=atoms, db_name=db_name)

# Run calculations for all structures that are not converged.
for row in db.select():
  atoms = row.toatoms()
  print("These are the atoms")
  print(atoms)


from clease import Evaluate
eva = Evaluate(settings=settings, scoring_scheme='loocv')

# scan different values of alpha and return the value of alpha that yields
# the lowest CV score
eva.set_fitting_scheme(fitting_scheme='l2')
#eva.set_fitting_scheme(fitting_scheme='l2', alpha=0.0001)
#alpha, cv = eva.alpha_CV(alpha_min=1E-6, alpha_max=10.0, num_alpha=50)
alpha = eva.plot_CV(alpha_min=1E-6, alpha_max=10.0, num_alpha=50)
eva.set_fitting_scheme(fitting_scheme='l2', alpha=alpha)
#print(eva.__dir__())
#print(eva.scheme)
#print(eva.cf_matr`ix)
#print(eva.e_dft)
eva.fit()

# plot ECI values
#eva.plot_ECI()
import clease.plot_post_process as pp
fig = pp.plot_fit(eva)
plt.show()

# plot ECI values
fig = pp.plot_eci(eva)
plt.show()
# save a dictionary containing cluster names and their ECIs
eva.save_eci(fname='eci_2023')
eci_file=open('eci_2023.json', 'r')
eci = json.load(eci_file)
print(eci)
from clease.calculator import attach_calculator
atoms = settings.atoms.copy()*(7, 7, 7)
atoms = attach_calculator(settings, atoms=atoms, eci=eci)


for i in range(0,len(atoms),1):
    if atoms[i].symbol=='Ti' and i%2 == 0:
        atoms[i].symbol='Zr'


print("These are the atoms")
print(atoms)
from clease.montecarlo import Montecarlo
from clease.montecarlo.observers import EnergyEvolution
from clease.montecarlo import Montecarlo
T = 2000
nsteps=1372000
mc = Montecarlo(atoms, T)

from clease.montecarlo.constraints import FixedElement
cnst = FixedElement('O')
mc.generator.add_constraint(cnst)

from clease.montecarlo.observers import EnergyEvolution, Snapshot, CorrelationFunctionObserver, AcceptanceRate
obs = EnergyEvolution(mc)
mc.attach(obs, interval=1000)
corr = CorrelationFunctionObserver(atoms.calc)
mc.attach(corr)
accrate = AcceptanceRate()
mc.attach(accrate)
snap = Snapshot(atoms, fname='snapshot')
mc.attach(snap, interval=137200)


mc.run(steps=nsteps)
energies = obs.energies
thermo = mc.get_thermodynamic_quantities()

f = open('MC_2000','w')
print(thermo, file=f)
print(energies, sep='\n', file=f)
f.close()

for i in range(2000, 0, -50):
    T = i
    #atoms = vasp.read_vasp(POSCARstring)
    #atoms = attach_calculator(settings, atoms=atoms, eci=eci)
    mc = Montecarlo(atoms, T)
    mc.generator.add_constraint(cnst)
    obs = EnergyEvolution(mc)
    mc.attach(obs, interval=1000)
    snap = Snapshot(atoms, fname='snapshot')
    mc.attach(snap, interval=137200)
    corr = CorrelationFunctionObserver(atoms.calc)
    mc.attach(corr)
    accrate = AcceptanceRate()
    mc.attach(accrate)
    
    mc.run(steps=137200)
    mc.run(steps=nsteps)

    energies = obs.energies
    rate=accrate.get_averages()
    thermo = mc.get_thermodynamic_quantities()
    corr_dict=corr.get_averages()
    
    from ase.build import sort
    from ase.io import vasp, db
    import pickle

    sorted_atoms=sort(atoms)
    vasp.write_vasp('POSCAR{}.vasp'.format(i), sorted_atoms, direct=False, wrap=False)
    POSCARstring='POSCAR{}.vasp'.format(i)
    
    f = open('MC_{}'.format(i),'w')
    print(thermo, file=f)
    print(energies, sep='\n', file=f)
    f.close()


