---
title: "Velocity Verlet Algorithm"
teaching: 0
exercises: 0
questions:
- "How do we practically solve Newton's equations of motion?"
objectives:
- "To Demonstrate the implementation of the Velocity-Verlet algorithm."
keypoints:
- "Velocity-Verlet algorithm provides a simple and stable numerical solution to Newton's equations of motion.  We can validate our implementation against the exactly-solvable dynamics of a classical harmonic oscillator."
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
    {: .language-python}
