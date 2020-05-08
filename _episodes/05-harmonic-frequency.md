---
title: "Computing the harmonic frequency"
teaching: 15
exercises: 15
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

### Harmonic Frequency
You might have learned that the Harmonic Oscillator potential, which is a reasonable model for the vibrational motion of diatomic atomcs near their equilibrium bond length, is given by

$$ V(r) = \frac{1}{2} k (r-r_{eq})^2 + V(r_{eq}) $$

and that the vibrational frequency of the molecule within the Harmonic oscillator model is given by

$$ \nu = \frac{1}{2\pi}\sqrt{\frac{k}{\mu}} $$

where $$ \mu $$ is the reduced mass of the molecule and $$ k $$ is known as the force constant, which we computed 
in the previous lesson.

The reduced mass of HF is defined as
$$ \mu = \frac{m_H \cdot m_F}{m_H + m_F}, $$

where $$ m_H $$ and $$ m_F $$ are the masses of Hydrogen and Fluoride, respectively.

Let's go ahead utilize the force constants from the last lesson to estimate the potential energy surface within the Harmonic approximation.  We can also estimate the vibrational frequency within the Harmonic approximation usin the force constant and the reduced mass.

{% include links.md %}
```
### define harmonic potential for each level of theory
RHF_Harm_Pot = RHF_k*(r_fine-RHF_Req)**2 + RHF_E_Spline(RHF_Req)
MP2_Harm_Pot = MP2_k*(r_fine-MP2_Req)**2 + MP2_E_Spline(MP2_Req)
CCSD_Harm_Pot = CCSD_k*(r_fine-CCSD_Req)**2 + CCSD_E_Spline(CCSD_Req)


### plot the Harmonic potential! 
plt.plot(r_fine, RHF_Harm_Pot, 'red')
plt.plot(r_fine, MP2_Harm_Pot, 'green')
plt.plot(r_fine, CCSD_Harm_Pot, 'blue')
plt.show()

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
