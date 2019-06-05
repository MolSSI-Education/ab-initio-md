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
