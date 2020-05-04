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

### Validating Velocity-Verlet algorithm with the Harmonic Oscillator
Newton's equation of motion can be solved analytically for the Harmonic oscillator, and we can use this fact to validate our Velocity-Verlet algorithm (which provides an approximate solution to Newton's equation of motion for arbitrary potentials). That is, the vibrational motion of a diatomic subject to a Harmonic potential predicted by the Velocity-Verlet algorithm should closely match the analytical solution. Analytically, the bond length as a function of time for a diatomic experiencing a harmonic potential is given by

$$ r(t) = A \: {\rm sin}\left(\sqrt{\frac{k}{\mu}} t + \phi \right) + r_{eq}, $$


where 

$$ A = \frac{r(0)}{{\rm sin}(\phi)} $$, 

$$ r(0) $$ is the initial separation, and $$ \phi $$ is the initial phase of the cycle; note that corresponding to this initial separation is an initial velocity given by

$$ v(0) = A \: \sqrt{\frac{k}{\mu}} {\rm cos}\left( \phi \right).  $$

Let's define a function harmonicposition that takes arguments of $$ \sqrt{\frac{k}{\mu}} $$ (om), $$ A $$ (amp), $$ \phi $$ (phase), $$ r_{eq} $$ (req), and time (t), and returns the separation.


