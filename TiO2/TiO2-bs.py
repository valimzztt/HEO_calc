from ase.calculators.espresso import Espresso
from ase.db import connect
from clease.tools import update_db

pseudo_dir = './pseudos/'

pseudopotentials = {'Ti': pseudo_dir + 'Ti.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                    'Zr': pseudo_dir + 'Zr.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                    'O': pseudo_dir + 'O.pbesol-n-kjpaw_psl.1.0.0.UPF'}

input_data = {
        'control': {
           'calculation': 'relax',
           'restart_mode': 'from_scratch',
            'outdir': './out',
           'prefix': 'tutorial',        
            },
        'system': {
           'ecutwfc': 40.424239707,
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
  print("This value did not converge")
  atoms = row.toatoms()
  atoms.calc = calc
  atoms.get_potential_energy()
  update_db(uid_initial=row.id, final_struct=atoms, db_name=db_name)
print("We have successfully calulated all the atoms energies")

for row in db.select(id=2):
  print("Iterating over database")
  atoms = row.toatoms()
  print("These are the atoms")
  print(atoms)
  atoms.calc = calc
  #energy = atoms.get_potential_energy()
  #print("The atomic energy is equal to")
  #print(energy)
