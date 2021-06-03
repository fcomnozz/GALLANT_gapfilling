# Plastic_degradation_gap_filling
Repository for the development of a gap filling module for Gallant. 
This module contains two differentiated submoduleds:
1) Homology gap filling: it takes a query SBML model and a list of sorted templates (SBML too), and performs gap filling based on COBRA algorithm. 
2) ModelSEED access: it reconstructs models gived a set of PATRIC IDs, retrieves them, performs modelSEED gap filling and translates them into BiGG nomenclature.
