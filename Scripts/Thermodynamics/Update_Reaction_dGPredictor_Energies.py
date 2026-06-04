#!/usr/bin/env python
import os,sys,json,glob
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
from Estimate_Reaction_Reversibility import reversibility_from_energy

# dGPredictor (Wang et al. 2021, Metab Eng) predicts reaction dG directly from a
# group decomposition + ML model, output in kJ/mol. Predictions are staged as
# raw JSON in Biochemistry/Thermodynamics/dGPredictor/json_files/, keyed
# ModelSEED-rxn -> KEGG-R-id -> {dG_mean, dG_uncer}.
#
# This script does two independent things:
#
#  1. ADDITIVE record (every predicted reaction): it stores the dGPredictor
#     estimate in the JSON thermodynamics dict under 'dGPredictor' as
#     [energy, error, operator] for EVERY reaction dGPredictor predicts,
#     sitting next to the Group-Contribution / eQuilibrator records rather than
#     replacing them. The operator is this estimate's own thermodynamic
#     direction.
#
#  2. CANONICAL gap-fill (only where needed): it additionally writes the
#     dGPredictor energy into the top-level deltag / deltagerr columns ONLY for
#     reactions that currently have NO Group-Contribution / eQuilibrator
#     estimate (deltag == sentinel 10000000), leaving the well-validated GC/eQ
#     canonical values untouched. Gap-filled reactions are tagged 'DGP' in notes
#     (the reaction-level provenance mechanism).

KJ_PER_KCAL = 4.184
SENTINEL = 10000000

label = "dGPredictor"
thermo_root = os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/dGPredictor/json_files/"

# rxn -> (kcal mean, kcal err) aggregated across its KEGG ids
dgp = dict()
for path in sorted(glob.glob(thermo_root+"reaction_*_dG.json")):
    with open(path) as fh:
        obj = json.load(fh)
    for rxn, kegg_map in obj.items():
        means=list(); uncs=list()
        if(isinstance(kegg_map, dict)):
            for kegg, payload in kegg_map.items():
                if(isinstance(payload, dict) and "dG_mean" in payload):
                    m = payload["dG_mean"]
                    u = payload.get("dG_uncer", 0.0)
                    if(isinstance(m,(int,float)) and m==m and abs(m)!=float("inf")):
                        means.append(m)
                        uncs.append(u if (isinstance(u,(int,float)) and u==u) else 0.0)
        if(means):
            dg_kcal = round(sum(means)/len(means)/KJ_PER_KCAL, 2)
            err_kcal = round(sum(uncs)/len(uncs)/KJ_PER_KCAL, 2)
            dgp[rxn] = (dg_kcal, err_kcal)

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

stored=0
filled=0
for rxn in sorted(reactions_dict.keys()):
    robj = reactions_dict[rxn]

    if(robj.get('status') == 'EMPTY'):
        continue

    if(rxn not in dgp):
        continue

    (dg_kcal, err_kcal) = dgp[rxn]

    # 1. ADDITIVE: record the dGPredictor estimate alongside any GC/eQ records,
    #    as [energy, error, operator]. Done for every predicted reaction.
    operator = reversibility_from_energy(robj, dg_kcal, err_kcal)
    if(not isinstance(robj.get('thermodynamics'), dict)):
        robj['thermodynamics'] = dict()
    robj['thermodynamics'][label] = [dg_kcal, err_kcal, operator]
    stored+=1

    # 2. CANONICAL gap-fill: only set the top-level deltag/deltagerr (and tag
    #    'DGP') when the reaction has NO GC/eQ estimate. GC/eQ canonical values
    #    are left untouched.
    dg = robj.get('deltag')
    has_value = isinstance(dg,(int,float)) and abs(dg) < SENTINEL
    if(has_value):
        continue

    robj['deltag'] = dg_kcal
    robj['deltagerr'] = err_kcal

    notes = robj.get('notes')
    if(not isinstance(notes, list)):
        notes = []
    if('DGP' not in notes):
        notes.append('DGP')
    robj['notes'] = notes

    filled+=1

print("dGPredictor reactions available: "+str(len(dgp)))
print("dGPredictor records stored additively: "+str(stored))
print("Gap-filled reactions (no prior GC/eQ deltag): "+str(filled))
print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
