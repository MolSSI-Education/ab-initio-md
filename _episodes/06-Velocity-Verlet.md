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

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

 <script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>

### Numerically solving Newton's equation of motion
If the acceleration, position, and velocity of the bond stretch coordinate are known at some instant in time $$ t_i $$, then the position and velocity can be estimated at some later time $$ t_{i+1} = t_i + \Delta t $$:

$$ r(t_i + \Delta t) = r(t_i) + v(t_i)\Delta t + \frac{1}{2}a(t_i)\Delta t^2 $$

and

$$ v(t_i + \Delta t) = v(t_i) + \frac{1}{2} \left(a(t_i) + a(t_i + \Delta t)  \right) \Delta t. $$

This prescription for updating the velocities and positions is known as the Velocity-Verlet algorithm.
Note that we need to perform 2 force evaluations per Velocity-Verlet iteration: one corresponding to position $$ r(t_i) $$ to update the position, and then a second time at the updated position $$ r(t_i + \Delta t) $$ to complete the velocity update.

We will create a function called Velocity_Verlet that takes the arguments r_curr, v_curr, mu, force_spline, and timestep and returns a 2-element array containing the updated position (r) and velocity (v) value.


{% include links.md %}

```
### Function that implements the Velocity-Verlet Algorithm
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
