---
title: "Introduction"
teaching: 30 min
exercises: 0
questions:
- "How can we simulate the vibrational motion of molecules subject to realistic intramolecular forces?"
objectives:
- "We can solve Newton's equations of motion for relative motion of atoms subject to intramolecular forces."
keypoints:
- "We can use quantum chemistry to obtain a potential energy surface, from which realistic intramolecular forces 
can be derived.  We can then perform molecular dynamics simulations to solve Newton's equation of motion for the relative motion of atoms subject to realistic intramolecular forces."
---
{% include links.md %}

Molecular Dynamics simulations provide a powerful tool for studying the motion of atomic, molecular, and nanoscale systems 
(including biomolecular systems).  In molecular dynamics simulation, the motion of the atoms are treated with classical mechanics; this means that the positions and velocities/momenta of the atoms are precisely specified throughout the simulation, and these quantities are updated using [Newton's laws of motion](https://en.wikipedia.org/wiki/Newton%27s_laws_of_motion) (or equivalent formulations).  

The motion of the atoms over the duration of a molecular dynamics simulation is known as a trajectory.  The analysis of 
molecular dynamics trajectories can provide information about thermodynamic and kinetic quantities relevant for your system, and can also yield valuable insights into the mechanism of processes which occur over time.  Thus molecular dynamics simulations are employed to study an enormous range of chemical, material, and biological systems.  

A brief overview of different applications of molecular dynamics simulations, as well as an introduction to the key working
equations can be found [here](https://github.com/FoleyLab/2019-Tapia-MolSSI/blob/master/Further_Reading/MD_MainPaper.pdf), with some additional details [here](https://github.com/FoleyLab/2019-Tapia-MolSSI/blob/master/Further_Reading/MD_TheoreticalBackground.pdf), as well as to references cited therein.

Forces come from _ab_ _initio_ potential energy surface

Potential energy surface is the ground-state electronic energy of a molecule at different geometries
