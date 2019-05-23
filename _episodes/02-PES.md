
---
title: "_ab_ _initio_ potential energy surface"
teaching: 0
exercises: 0
questions:
- "How can we simulate the vibrational motion of a molecule?"
objectives:
- "We can solve Newton's equations of motion for relative motion of atoms subject to intramolecular forces."
keypoints:
- "We can use quantum chemistry to obtain accurate intramolecular forces, and we can use molecular dynamics simulations 
to solve Newton's equations of motion to study the vibrational motion moving according to intramolecular force."
---
FIXME

{% include links.md %}

``` {: .language-python}
import numpy as np
import psi4
from matplotlib import pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
```

``` {: .language-python}
#This is an example of the code we would expect a student to write.  
#If they use the exact same values for rvals, they don't find the well of the PES.  
#The key is that they realize they need to use more R values to get the right shape graph.
mol_tmpl = """H
F 1 **R**"""
r_array = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3]
molecules =[]
HF_E_array = []
MP2_E_array = []
CCSD_E_array = []
for r in r_array:
    molecule = psi4.geometry(mol_tmpl.replace("**R**", str(r)))
    molecules.append(molecule)
for mol in molecules:
    energy = psi4.energy("SCF/cc-pVDZ", molecule=mol)
    HF_E_array.append(energy)
    energy = psi4.energy("MP2/cc-pVDZ", molecule=mol)
    MP2_E_array.append(energy)
    energy = psi4.energy("CCSD/cc-pVDZ",molecule=mol)
    CCSD_E_array.append(energy)
    
plt.plot(r_array,HF_E_array,'*')
plt.plot(r_array,MP2_E_array,'*')
plt.plot(r_array,CCSD_E_array,'*')
```
