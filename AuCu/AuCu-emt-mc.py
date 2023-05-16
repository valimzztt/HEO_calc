import numpy as np

# 1) Specify the concentration ranges of species
from clease.settings import Concentration
conc = Concentration(basis_elements=[['Au', 'Cu']], A_lb=[[2, 0]], b_lb=[1])
# 2) Stochiometric Constraints: Binary System With One Basis
basis_elements = [['Au', 'Cu']]
A_eq = [[1.0, -1.0]]
b_eq = [0.0]
conc = Concentration(basis_elements=basis_elements, A_eq=A_eq, b_eq=b_eq)
for i in range(10):
    x = conc.get_random_concentration([20])
    assert np.abs(x[0] - x[1]) < 1e-10
# 3) Specify CE settings
from clease.settings import CEBulk
settings = CEBulk(crystalstructure='fcc',
                  a=3.8,
                  supercell_factor=64,
                  concentration=conc,
                  db_name="aucu.db",
                  max_cluster_dia=[6.0, 4.5, 4.5])


# 4) Generating new structures
from ase.db import connect
from clease.structgen import NewStructures
ns = NewStructures(settings, generation_number=0,
                   struct_per_gen=10)

# Get template with the cell size = 3x3x3
atoms = connect('aucu.db').get(id=10).toatoms()
ns.generate_random_structures(atoms)

# 5) Running calculations on generated structures
from ase.calculators.emt import EMT
calc = EMT()
from ase.db import connect
from clease.tools import update_db
db_name = "aucu.db"
db = connect(db_name)
# (Run calculations for all structures that are not converged.)
for row in db.select(converged=False):
  atoms = row.toatoms()
  atoms.calc = calc
  atoms.get_potential_energy()
  update_db(uid_initial=row.id, final_struct=atoms, db_name=db_name)

# 6) Evaluating the CE model
from clease import Evaluate
import clease.plot_post_process as pp
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
eva = Evaluate(settings=settings, scoring_scheme='k-fold', nsplits=10)
# scan different values of alpha and return the value of alpha that yields
# the lowest CV score
eva.set_fitting_scheme(fitting_scheme='l1')
alpha = eva.plot_CV(alpha_min=1E-7, alpha_max=1.0, num_alpha=50)
# set the alpha value with the one found above, and fit data using it.
print("Optimal value of alpha is equal to")
print(alpha)
eva.set_fitting_scheme(fitting_scheme='l1', alpha=alpha)
eva.fit()  
# Run the fit with these settings.
fig = pp.plot_fit(eva)
plt.show()
# plot ECI values (J_alpha)
fig = pp.plot_eci(eva)
plt.show()
# save a dictionary containing cluster names and their ECIs
eva.save_eci(fname='eci_l1')



# 7) Generating structures for further training
from clease import NewStructures
ns = NewStructures(settings, generation_number=1, struct_per_gen=10)
ns.generate_probe_structure()

#After these new structures are created, rerun part number 5 of main.py script

#Next step to follow after having constructed the CE model
# 1) Monte Carlo Sampling
# CLEASE currently support two ensembles for Monte Carlo sampling: canonical and semi-grand canonical ensembles.

from clease.settings import CEBulk, Concentration
conc = Concentration(basis_elements=[['Au', 'Cu']])
settings = CEBulk(crystalstructure='fcc',
                  a=3.8,
                  supercell_factor=27,
                  concentration=conc,
                  db_name="aucu.db",
                  max_cluster_dia=[6.0, 5.0])
# specify as set of ECIs
eci = {'c0': -1.0, 'c1_0': 0.1, 'c2_d0000_0_00': -0.2}
from clease.calculator import attach_calculator
atoms = settings.atoms.copy()*(5, 5, 5)
atoms = attach_calculator(settings, atoms=atoms, eci=eci)
print(atoms.ase_objtype)
atoms[0].symbol = 'Cu'
atoms[1].symbol = 'Cu'
atoms[2].symbol = 'Cu'

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





