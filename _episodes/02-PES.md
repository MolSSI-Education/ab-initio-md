
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

``` 
import numpy as np
import psi4
from matplotlib import pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
```
{: .language-python}

``` 
### template for the z-matrix
mol_tmpl = """H
F 1 **R**"""
### array of bond-lengths in anstromgs
r_array = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3]
### array for different instances of the HF molecule
molecules =[]
### array for the different RHF energies for different HF bond-lengths
HF_E_array = []
### array for the different MP2 energies for different HF bond-lengths
MP2_E_array = []
### array for the different CCSD energies for different HF bond-lengths
CCSD_E_array = []

### loop over the different bond-lengths, create different instances
### of HF molecule
for r in r_array:
    molecule = psi4.geometry(mol_tmpl.replace("**R**", str(r)))
    molecules.append(molecule)
    
### loop over instances of molecules, compute the RHF, MP2, and CCSD
### energies and store them in their respective arrays
for mol in molecules:
    energy = psi4.energy("SCF/cc-pVDZ", molecule=mol)
    HF_E_array.append(energy)
    energy = psi4.energy("MP2/cc-pVDZ", molecule=mol)
    MP2_E_array.append(energy)
    energy = psi4.energy("CCSD/cc-pVDZ",molecule=mol)
    CCSD_E_array.append(energy)

### Plot the 3 different PES
plt.plot(r_array,HF_E_array,'r*', label='RHF')
plt.plot(r_array,MP2_E_array,'b*', label='MP2')
plt.plot(r_array,CCSD_E_array,'g*', label='CCSD')
```
{: .language-python}

We will now define two arrays: r_array will be an array of values for the $HF$ bond length and E_array will hold the electronic energy values corresponding to each separation.  

``` 
### use cubic spline interpolation
order = 3
### form the interpolator object for the data
for i in range(0,len(r_array)):
    val = r_array[i]/0.529
    r_array[i] = val
    
sE = InterpolatedUnivariateSpline(r_array, E_array, k=order)
### form a much finer grid
r_fine = np.linspace(0.5/0.529,2.3/0.529,200)
### compute the interpolated/extrapolated values for E on this grid
E_fine = sE(r_fine)
### plot the interpolated data
plt.plot(r_fine, E_fine, 'blue')
plt.show()

### take the derivative of potential
fE = sE.derivative()
### force is the negative of the derivative
F_fine = -1*fE(r_fine)

### plot the forces
plt.plot(r_fine, F_fine, 'black')
plt.show()
```
{: .language-python}
