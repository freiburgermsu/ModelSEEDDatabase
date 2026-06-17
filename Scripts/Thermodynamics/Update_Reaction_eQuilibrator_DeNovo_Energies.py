#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
from Estimate_Reaction_Reversibility import reversibility_from_energy

# Writes the de-novo eQuilibrator estimates (from
# Retrieve_eQuilibrator_DeNovo_Reaction_Energies.py) ADDITIVELY into each
# reaction's thermodynamics dict under a DISTINCT source key,
# 'eQuilibrator-deNovo', as [energy, error, operator]. It is kept separate from
# the cache-derived 'eQuilibrator' source on purpose: de-novo estimates are
# lower-confidence (group-contribution on novel structures; default protonation
# when ChemAxon is unavailable) and must never be conflated or double-counted
# with the cache values. The canonical top-level deltag/deltagerr/reversibility
# are NOT modified.

label = "eQuilibrator-deNovo"
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

thermo_root = os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/"
file_name = thermo_root+'eQuilibrator/MetaNetX_DeNovo_Reaction_Energies.tbl'

eq = dict()
if(os.path.exists(file_name)):
    with open(file_name) as fh:
        for line in fh:
            a = line.rstrip("\n").split('\t')
            if(len(a) < 3 or 'energy' in a[1] or a[1] == 'nan'):
                continue
            eq[a[0]] = [float("{0:.2f}".format(float(a[1]))), float("{0:.2f}".format(float(a[2])))]

stored = 0
for rxn in sorted(reactions_dict.keys()):
    if(rxn not in eq):
        continue
    robj = reactions_dict[rxn]
    if(robj.get('status') == 'EMPTY'):
        continue

    (dg_val, dge_val) = eq[rxn]
    operator = reversibility_from_energy(robj, dg_val, dge_val)
    if(not isinstance(robj.get('thermodynamics'), dict)):
        robj['thermodynamics'] = dict()
    robj['thermodynamics'][label] = [dg_val, dge_val, operator]
    stored += 1

print("de-novo eQuilibrator reaction records: "+str(len(eq)))
print("stored additively under '"+label+"': "+str(stored))
print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
