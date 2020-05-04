---
title: "_ab_ _initio_ potential energy surface"
teaching: 15 minutes
exercises: 15 minutes
questions:
- "How do we obtain realistic intramolecular forces?"
objectives:
- "To provide a simple demonstration of the computation of potential energy surfaces for a diatomic molecule and subsequent determination of intramolecular forces."
keypoints:
- "Psi4numpy provides an easy framework for computing and manipulating quantum chemical properties, including molecular energies required for potential energy surfaces.  Scipy's interpolation libraries can facilitate the fitting of potential
energy surfaces, and for automatically computing derivatives of these surfaces to obtain intramolecular forces."
---

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

Now that you have the raw data, we will interpolate this data using cubic splines.  This will permit us to 
estimate the potential energy at any arbitrary separation between 0.5 and 2.3 Angstroms (roughly 
1 and 4.3 a.u.) with fairly high confidence for each level of theory, and will also allow us to estimate the corresponding
forces 
(the negative of the derivative of the PES with respect to separation)
at any separation between 1.0 and 4.3 a.u. since the derivative of cubic splines are readily available.

First, we wwill interpolate the energies for each level of theory:

``` 
### use cubic spline interpolation
order = 3

### get separation vector in atomic units
r_array_au = 1.89*r_array 

### spline for RHF Energy
RHF_E_Spline = InterpolatedUnivariateSpline(r_array_au, HF_E_array, k=order)

### spline for MP2 Energy
MP2_E_Spline = InterpolatedUnivariateSpline(r_array_au, MP2_E_array, k=order)

### spline for CCSD Energy
CCSD_E_Spline = InterpolatedUnivariateSpline(r_array_au, CCSD_E_array, k=order)


### form a much finer grid
r_fine = np.linspace(0.5/0.529,2.3/0.529,200)

### compute the interpolated/extrapolated values for RHF Energy on this grid
RHF_E_fine = RHF_E_Spline(r_fine)

### compute the interpolated/extrapolated values for RHF Energy on this grid
MP2_E_fine = MP2_E_Spline(r_fine)

### compute the interpolated/extrapolated values for RHF Energy on this grid
CCSD_E_fine = CCSD_E_Spline(r_fine)


### plot the interpolated data to check our fit!
plt.plot(r_fine, RHF_E_fine, 'red', r_array_au, HF_E_array, 'r*')
plt.plot(r_fine, MP2_E_fine, 'green', r_array_au, MP2_E_array, 'g*')
plt.plot(r_fine, CCSD_E_fine, 'blue', r_array_au, CCSD_E_array, 'b*')

plt.show()
```
{: .language-python}

Next, we can compute the forces!  Recall that the force comes from the negative derivative of the potential,

$$ F(r) = -\frac{d}{dr} V(r). $$

```
### take the derivative of the potential to get the negative of the force from RHF
RHF_Force = RHF_E_Spline.derivative() 

### negative of the force from MP2
MP2_Force = MP2_E_Spline.derivative()

### negative of the force from CCSD
CCSD_Force = CCSD_E_Spline.derivative()

### let's plot the forces for each level of theory!

### plot the forces... note we need to multiply by -1 since the spline
### derivative gave us the negative of the force!
plt.plot(r_fine, -1*RHF_Force(r_fine), 'red')
plt.plot(r_fine, -1*MP2_Force(r_fine), 'green')
plt.plot(r_fine, -1*CCSD_Force(r_fine), 'blue')
plt.show()
```
{: .language-python}
