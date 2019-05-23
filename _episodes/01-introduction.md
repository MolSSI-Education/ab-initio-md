---
title: "Introduction"
teaching: 0
exercises: 0
questions:
- "How can we simulate the vibrational motion of molecules subject to realistic intramolecular forces?"
objectives:
- "We can solve Newton's equations of motion for relative motion of atoms subject to intramolecular forces."
keypoints:
- "We can use quantum chemistry to obtain a potential energy surface, from which realistic intramolecular forces 
can be derived.  We can then perform molecular dynamics simulations to solve Newton's equation of motion for the relative motion of atoms subject to realistic intramolecular forces."
---
FIXME

{% include links.md %}

Now that you have the raw data, we will interpolate this data using cubic splines.  This will permit us to 
estimate the potential energy at any arbitrary separation between 0.5 and 3.5 Angstroms (roughly 
1 and 5.8 a.u.) with fairly high confidence, and will also allow us to estimate the force 
\begin{equation}
F(r) = -\frac{d}{dr} V(r)
\end{equation}
at any separation between 1.0 and 5.8 a.u. since the derivative of cubic splines are readily available.
