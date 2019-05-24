---
title: "_ab_ _initio_ potential energy surface"
teaching: 0
exercises: 0
questions:
- "How do we obtain realistic intramolecular forces?"
objectives:
- "To provide a simple demonstration of the computation of potential energy surfaces for a diatomic molecule and subsequent determination of intramolecular forces."
keypoints:
- "Psi4numpy provides an easy framework for computing and manipulating quantum chemical properties, including molecular energies required for potential energy surfaces.  Scipy's interpolation libraries can facilitate the fitting of potential
energy surfaces, and for automatically computing derivatives of these surfaces to obtain intramolecular forces."
---

{% include links.md %}

```
def Velocity_Verlet(r_curr, v_curr, mu, f_interp, dt):
    ### get acceleration at current time
    a_curr = -1*f_interp(r_curr)/mu
    
    ### use current acceleration and velocity to update position
    r_fut = r_curr + v_curr * dt + 0.5 * a_curr * dt**2
    
    ### use r_fut to get future acceleration a_fut
    a_fut = -1*f_interp(r_fut)/mu
    ### use current and future acceleration to get future velocity v_fut
    v_fut = v_curr + 0.5*(a_curr + a_fut) * dt
    
    result = [r_fut, v_fut]
    
    return result
    ```
    
