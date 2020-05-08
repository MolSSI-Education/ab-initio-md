---
title: "Introduction"
teaching: 30 min
exercises: 0
questions:
- "How can we simulate the motion of molecules subject to realistic intramolecular forces?"
objectives:
- "We can solve Newton's equations of motion for relative motion of atoms subject to intramolecular forces."
keypoints:
- "We can use quantum chemistry to obtain a potential energy surface, from which realistic intramolecular forces 
can be derived.  We can then perform molecular dynamics simulations to solve Newton's equation of motion for the relative motion of atoms subject to realistic intramolecular forces."
---
{% include links.md %}

<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML"> </script> <script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>

Molecular Dynamics simulations provide a powerful tool for studying the motion of atomic, molecular, and nanoscale systems 
(including biomolecular systems).  In molecular dynamics simulation, the motion of the atoms are treated with classical mechanics; this means that the positions and velocities/momenta of the atoms are precisely specified throughout the simulation, and these quantities are updated using [Newton's laws of motion](https://en.wikipedia.org/wiki/Newton%27s_laws_of_motion) (or equivalent formulations).  

A brief overview of different applications of molecular dynamics simulations, as well as an introduction to the key working
equations can be found [here](https://github.com/FoleyLab/2019-Tapia-MolSSI/blob/master/Further_Reading/MD_MainPaper.pdf), with some additional details [here](https://github.com/FoleyLab/2019-Tapia-MolSSI/blob/master/Further_Reading/MD_TheoreticalBackground.pdf), as well as to references cited therein.  An even more concise summary of this excercise can be found in [these slides](https://github.com/FoleyLab/2019-Tapia-MolSSI/blob/master/Overview.pdf).

The motion of the atoms over the duration of a molecular dynamics simulation is known as a trajectory.  The analysis of 
molecular dynamics trajectories can provide information about thermodynamic and kinetic quantities relevant for your system, and can also yield valuable insights into the mechanism of processes which occur over time.  Thus molecular dynamics simulations are employed to study an enormous range of chemical, material, and biological systems.  

One of the central quantities in a molecular dynamics simulation comes from the inter-particle forces.  From Newton's second law, we know that the acceleration a particle feels is related to the force it experiences;

$$ \vec{F} = m \vec{a}. $$

The acceleration of each particle in turn determines how the position and momentum of each particle changes.  Therefore, the trajectories one observes in a molecular dynamics simulation is governed by the forces experienced by the particles being simulated.

In this exercise, we will compute the forces experienced by our particles using the tools of quantum chemistry.  The resulting molecular dynamics simulation is often referred to as an _ab_ _initio_ molecular dynamics simulation, since the underlying forces are determined from _first_ _principles_, i.e. through quantum mechanics.  To the extent that our _ab_ _initio_ method is accurate, we will obtain accurate forces to recover realistic dynamics of the system we are studying.   


