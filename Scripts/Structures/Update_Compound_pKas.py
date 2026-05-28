#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["SMILE","InChI","InChIKey"],["ModelSEED"])
aliases_dict = compounds_helper.loadMSAliases()

# Load pKas and pKbs from the post-A1 layout: <db>/pkas/<tool>_<ver>.tsv.
# Iteration order is the priority cascade KEGG > MetaCyc > ChEBI > Rhea.
# First DB with a hit on any of the compound's aliases wins.
#
# ChEBI and Rhea were previously unused because of ID-format mismatches
# (ChEBI: 'CHEBI_15377' in pKa file vs '15377' in aliases; Rhea:
# 'POLYMER_X' vs 'POLYMER:X'). Compounds.loadPerSourcePkas normalizes
# both at ingestion now, see sources.yaml for the migration history.
# Adding ChEBI to the cascade unlocks ~3,800 additional pKa
# attributions; Rhea's polymer pKa data is fully shadowed by primary
# sources so its inclusion is for completeness rather than coverage.
PKA_DBS = ["KEGG", "MetaCyc", "ChEBI", "Rhea"]
per_source_pkas = compounds_helper.loadPerSourcePkas(PKA_DBS)
cpd_pKab_dict = dict()
for (db, ext_id), entry in per_source_pkas.items():
    if(ext_id not in cpd_pKab_dict):
        cpd_pKab_dict[ext_id] = dict()
    for kind, value in entry.items():
        cpd_pKab_dict[ext_id][kind] = value

# OPAM2 / MolGpKa pKas are a ModelSEED-compound-level source (keyed by cpd id,
# not by per-source external id) and take precedence over Marvin: where OPAM2
# has a prediction it OVERRIDES the Marvin value applied below; Marvin is
# retained for compounds OPAM2 does not cover. OPAM2 predicts pKa values only
# (no protonation/structure), so formula/charge are untouched. Built from the
# OPAM2 benchmark (Marvin-matched atom subset) by Scripts/Updates helpers.
OPAM2_PKA_FILE = os.path.dirname(__file__)+"/../../Biochemistry/Structures/ModelSEED/pkas/opam2_molgpka.tsv"
opam2_pkas = dict()  # cpd -> {'pKa':str, 'pKb':str}
if(os.path.exists(OPAM2_PKA_FILE)):
    with open(OPAM2_PKA_FILE) as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if(len(cols) < 3):
                continue
            opam2_pkas.setdefault(cols[0], dict())[cols[1]] = cols[2]

# We're removing all pKa and pKb before loading new ones
for cpd in compounds_dict:
    compounds_dict[cpd]['pka']=""
    compounds_dict[cpd]['pkb']=""

# We're only loading pKa/pKb for compounds that have an accepted unique structure in ModelSEED
for cpd in structures_dict:
    found=False
    for DB in PKA_DBS:
        if(found is True or DB not in aliases_dict[cpd]):
            continue

        for alias in aliases_dict[cpd][DB]:
            if(alias in cpd_pKab_dict):
                print(cpd,alias,cpd_pKab_dict[alias])
                if('pKa' in cpd_pKab_dict[alias]):
                    print(cpd,alias)
                    compounds_dict[cpd]['pka']=cpd_pKab_dict[alias]['pKa']
                else:
                    compounds_dict[cpd]['pka']=""

                if('pKb' in cpd_pKab_dict[alias]):
                    compounds_dict[cpd]['pkb']=cpd_pKab_dict[alias]['pKb']
                else:
                    compounds_dict[cpd]['pkb']=""

                # All structures for the same compound should be the same,
                # and so only need to process once
                found=True
                break

# OPAM2 override (primary pKa source). Applied only to compounds with an
# accepted unique structure, matching the Marvin gating above. A compound for
# which OPAM2 found only acidic (or only basic) atoms gets the complementary
# field cleared, reflecting the OPAM2 prediction rather than mixing tools.
opam2_applied=0
for cpd in structures_dict:
    if(cpd in opam2_pkas and cpd in compounds_dict):
        compounds_dict[cpd]['pka']=opam2_pkas[cpd].get('pKa',"")
        compounds_dict[cpd]['pkb']=opam2_pkas[cpd].get('pKb',"")
        opam2_applied+=1
print("Applied OPAM2 pKa overrides to "+str(opam2_applied)+" compounds")

print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
