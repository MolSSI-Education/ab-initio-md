---
title: "Choosing intial bond-length and velocity"
teaching: 0
exercises: 0
questions:
- "How do we initialize the position and velocity of the atoms for our molecular dynamics simulation?"
objectives:
- "To demonstrate a reasonable strategy for initializing a molecular dynamics simulation of vibrational motion."
keypoints:
- "We can utilize the virial theorem from quantum mechanics to provide a reasonable value of the classical velocity associated with vibrational motion."
---

{% include links.md %}
```
### define "ground-state" velocity for each level of theory
v_RHF = np.sqrt( np.sqrt(RHF_k/mu)/(2*mu))
v_MP2 = np.sqrt( np.sqrt(MP2_k/mu)/(2*mu))
v_CCSD = np.sqrt( np.sqrt(CCSD_k/mu)/(2*mu))


### get random position and velocity for RHF HF within a reasonable range
r_init = np.random.uniform(0.75*RHF_Req,2*RHF_Req)
v_init = np.random.uniform(-2*v_RHF,2*v_RHF)

### print initial position and velocity
print("Initial separation is ",r_init, "atomic units")
print("Initial velocity is   ",v_init, "atomic units")


### get initial force on the particle based on its separation
RHF_F_init = -1*RHF_Force(r_init)
print("Initial Force is ", RHF_F_init, "atomic units")
```
