---
title: "Force Constant for Harmonic Potential"
teaching: 15 minutes
exercises: 15 minutes
questions:
- "How do we define the potential within the Harmonic approxmation?"
objectives:
- "The second derivative at the equlibrium bond length is used to define the force constant, which defines the harmonic potential and force."
keypoints:
- "Each level of theory has a unique potential energy curve and might be expected to give a unique force constant."
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

 <script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>
{% include links.md %}

```
### get second derivative of potential energy curve... recall that we fit a spline to
### to the first derivative already and called that spline function X_Force, where
### X is either RHF, MP2, or CCSD

RHF_Curvature = RHF_Force.derivative()
MP2_Curvature = MP2_Force.derivative()
CCSD_Curvature = CCSD_Force.derivative()

### evaluate the second derivative at r_eq to get k
RHF_k = RHF_Curvature(RHF_Req)
MP2_k = MP2_Curvature(MP2_Req)
CCSD_k = CCSD_Curvature(CCSD_Req)

### Print force constants for each level of theory!
print("Hartree-Fock force constant is ",RHF_k," atomic units")
print("MP2 force constant is ",MP2_k," atomic units")
print("CCSD force constant is ",CCSD_k," atomic units")


### define harmonic potential for each level of theory
RHF_Harm_Pot = RHF_k*(r_fine-RHF_Req)**2 + RHF_E_Spline(RHF_Req)
MP2_Harm_Pot = MP2_k*(r_fine-MP2_Req)**2 + MP2_E_Spline(MP2_Req)
CCSD_Harm_Pot = CCSD_k*(r_fine-CCSD_Req)**2 + CCSD_E_Spline(CCSD_Req)


### plot!
plt.plot(r_fine, RHF_Harm_Pot, 'red')
plt.plot(r_fine, MP2_Harm_Pot, 'green')
plt.plot(r_fine, CCSD_Harm_Pot, 'blue')
plt.show()
```
{: .language-python}

## The "force field"

As stated above, all computational simulations depend on a mathematical model to describe the molecule. The full expression of this model is often called a **force field**  and is a function which describes the potential energy of a system. This force field describes the energy associated with molecular movements such as bond stretching, angle bending, or dihedral angle rotation. If you are a chemistry student, you have likely discussed many of these molecular motions in your classes (who can forget talking about the cis-trans isomerism of butane in organic chemistry class?). The force field describes the energies associated with these movements mathematically. The potential energy is commonly represented by the letter $$U$$, and the total potential energy or "force field" is the sum of terms related to bond stretching, angle bending, torsional rotation (that cis-trans isomerism), electrostatic interaction, and nonbonded Van der Waals interactions.

$$ U = U_{bond} + U_{angle} + U_{torsion} + U_{elec} + U_{vdw} $$
