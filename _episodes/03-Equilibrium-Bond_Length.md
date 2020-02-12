---
title: "Finding the equilibrium bond length"
teaching: 5 minutes
exercises: 10 minutes
questions:
- "How do we find the equilibrium bond length of our diatomic molecule?"
objectives:
- "The equilibrium bond length for a diatomic molecule is defined as the separation that minimizes the total electronic energy; we can identify this as the minimum on our potential energy curves."
keypoints:
- "Each level of theory has a unique potential energy curve and might be expected to give a unique equilibrium bond length.  We have fit each potential energy curve with splines, so there may be limitations to the precision in our equilibrium bond length estimate coming from the resultion of the points we fit our spline to, or from the fit of the spline itself; which do you think is more significant?"
---

{% include links.md %}
```
### Find Equilibrium Bond-Lengths for each level of theory
RHF_Req_idx = np.argmin(RHF_E_fine)
MP2_Req_idx = np.argmin(MP2_E_fine)
CCSD_Req_idx = np.argmin(CCSD_E_fine)

### find the value of the separation corresponding to that index
RHF_Req = r_fine[RHF_Req_idx]
MP2_Req = r_fine[MP2_Req_idx]
CCSD_Req = r_fine[CCSD_Req_idx]

### print equilibrium bond-lengths at each level of theory!
print(" Equilibrium bond length at RHF/cc-pVDZ level is ",RHF_Req)
print(" Equilibrium bond length at MP2/cc-pVDZ level is ",MP2_Req)
print(" Equilibrium bond lengthat CCSD/cc-pVDZ level is ",CCSD_Req)
```
{: .language-python}
