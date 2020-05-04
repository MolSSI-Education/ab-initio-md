---
title: "Computing the harmonic frequency"
teaching: 0
exercises: 0
questions:
- "How do we use information about the _ab_ _initio_ potential energy surface to estimate the vibrational frequency of a diatomic molecule??"
objectives:
- "We can use the harmonic approximation to estimate the vibrational frequency."
keypoints:
- "Within the harmonic oscillator approximation, the frequency is related to the force constant divided by the reduced mass.  The force constant is defined as the second derivative of the PES at the equilibrium bond length."
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

 <script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>

###. Harmonic Frequency
You might have learned that the Harmonic Oscillator potential, which is a reasonable model for the vibrational motion of diatomic atomcs near their equilibrium bond length, is given by

$$ V(r) = \frac{1}{2} k (r-r_{eq})^2 + V(r_{eq}) $$

and that the vibrational frequency of the molecule within the Harmonic oscillator model is given by

$$ \nu = \frac{1}{2\pi}\sqrt{\frac{k}{\mu}} $$

where $$ \mu $$ is the reduced mass of the molecule and $$ k $$ is known as the force constant, which as we saw 
in the previous lesson can be defined as

$$ k = \frac{d^2}{dr^2} V(r_{eq}).$$

The reduced mass of HF is defined as
$$ \mu = \frac{m_H \cdot m_F}{m_H + m_F}, $$

where $$ m_H $$ and $$ m_F $$ are the masses of Hydrogen and Fluoride, respectively.

Let's go ahead and get the force constants at each level of theory, print the values, and estimate the potential energy within the Harmonic approximation! Just like we were able to differentiate our PES splines to get a force spline, we can differentiate a force splines to get curvature splines (which we can call RHF_Curvature, MP2_Curvature, and CCSD_Curvature); the force constant will then be the curvature evaluated at the equlibrium bond length.

{% include links.md %}
```

### define reduced mass of HF as m_H * m_H /(m_F + m_H) where mass is in atomic units (electron mass = 1)
m_F = 34883.
m_H = 1836.
mu = (m_F * m_H)/(m_F + m_H)

### compute the fundamental frequency at each level of theory
RHF_nu = 1/(np.pi*2) * np.sqrt(RHF_k/mu)
MP2_nu = 1/(np.pi*2) * np.sqrt(MP2_k/mu)
CCSD_nu = 1/(np.pi*2) * np.sqrt(CCSD_k/mu)

### print the values in atomic units!
print("Vibrational frequency of HF at the RHF/cc-pVDZ level is ",RHF_nu," atomic units")
print("Vibrational frequency of HF at the MP2/cc-pVDZ level is ",MP2_nu," atomic units")
print("Vibrational frequency of HF at the CCSD/cc-pVDZ level is ",CCSD_nu," atomic units")
```
{: .language-python}
