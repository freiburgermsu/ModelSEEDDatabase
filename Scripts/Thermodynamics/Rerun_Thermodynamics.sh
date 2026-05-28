#!/bin/bash
./Update_Compound_GroupContribution_Energies.py
./Update_Reaction_GroupContribution_Energies.py
./Estimate_Reaction_Reversibility.py GC
./Update_Compound_eQuilibrator_Energies.py
./Update_Reaction_eQuilibrator_Energies.py
./Estimate_Reaction_Reversibility.py EQ
# Gap-fill reactions with no GC/eQ estimate using staged dGPredictor predictions
./Update_Reaction_dGPredictor_Energies.py
