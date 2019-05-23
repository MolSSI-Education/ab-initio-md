
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

# Computational Exercise 2-4: A Surface with so much potential
We are going to construct what is often referred to as an $ab$ $initio$ potential energy surface of the diatomic
molecule carbon monoxide.  That is, we are going to use various electronic structure theories (Hartree-Fock theory, Configuration Interaction theory, and Density Functional theory) to compute the electronic energy at different geometries of a simple diatomic molecule.  We will use the interpolation capabilities of scipy to simplify the evalution of the potential energy at separations for which we did not explicitly evaluate the electronic energy.  We will also use scipy to differentiate the interpolated potential energy surface to obtain the forces acting on the atoms at different separations.  

We will start by importing the necessary libraries:

`import numpy as np`

`import psi4`

`from matplotlib import pyplot as plt`

`from scipy.interpolate import InterpolatedUnivariateSpline`
