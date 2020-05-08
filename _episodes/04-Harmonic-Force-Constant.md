---
title: "Force Constant for Harmonic Potential"
teaching: 5 minutes
exercises: 10 minutes
questions:
- "How do we define the force constant of a diatomic molecule?"
objectives:
- "The second derivative at the equlibrium bond length is used to define the force constant."
keypoints:
- "Each level of theory has a unique potential energy curve and might be expected to give a unique force constant."
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

 <script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>
{% include links.md %}

## The force constant

The force constant is defined as the second derivative of the potential energy surface evaluated at the minimum,

$$ k = \frac{d^2}{dr^2} V(r_{eq}). $$

With the equilibrium bond-length and the force spline in hand, we can utilize the derivative method associated with our force splines and evaluate it at the equilibrium bond length to compute the force constants.

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

```
{: .language-python}


