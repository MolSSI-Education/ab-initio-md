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
Newton's equation of motion can be solved analytically for the Harmonic oscillator, and we can use this fact to validate our Velocity-Verlet algorithm. That is, the vibrational motion of a diatomic subject to a Harmonic potential predicted by the Velocity-Verlet algorithm should closely match the analytical solution. Analytically, the bond length as a function of time for a diatomic experiencing a harmonic potential is given by

$$ r(t) = A \: {\rm sin}\left( \sqrt{\frac{k}{\mu}} t + \phi \right) + r_{eq} $$ where 

$$ A = \frac{r(0)}{ {\rm sin}(\phi) } $$, 

$$ r(0) $$ is the initial separation, and $$ \phi $$ is the initial phase of the cycle; note that corresponding to this initial separation is an initial velocity given by

$$ v(0) = A \: \sqrt{\frac{k}{\mu}} {\rm cos}\left( \phi \right).  $$

Let's define a function harmonicposition that takes arguments of $$ \sqrt{\frac{k}{\mu}} $$ (om), $$ A $$ (amp), $$ \phi $$ (phase), $$ r_{eq} $$ (req), and time (t), and returns the separation.


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
    
''' Students will write this! '''
def harmonic_position(om, Amp, phase, req, time):   
    return  Amp * np.sin( om * time + phase ) + req
    

''' This will be pre-written! '''
### how many updates do you want to perform?
N_updates = 10000

### establish time-step for integration to be 0.02 atomic units... this is about 0.0005 femtoseconds
### so total time is 200000*0.02 atomic units of time which is ~9.6e-13 s, or 960 fs
dt = 0.1

### results from VV algorithm
hr_vs_t = np.zeros(N_updates)
hv_vs_t = np.zeros(N_updates)
### analytic result for r(t)
ar_vs_t = np.zeros(N_updates)
### array to store time in atomic units
t_array = np.zeros(N_updates)

### establish some constants relevant for analytic solution
### harmonic freq
om = np.sqrt(RHF_k/mu)
### initial displacement 
x0 = 0.2
### amplitude for analytic solution
Amp = x0/(np.sin(np.pi/4))
### initial velocity
v0 = Amp * om * np.cos(np.pi/4)

hr_vs_t[0] = RHF_Req+x0
hv_vs_t[0] = v0

### We need a spline object for the harmonic force to pass to the Velocity Verlet algorithm,
### let's get that now!
### spline for Harmonic potential using RHF_k
RHF_Harm_Pot_Spline = InterpolatedUnivariateSpline(r_fine, RHF_Harm_Pot, k=3)
### RHF harmonic force
RHF_Harm_Force = RHF_Harm_Pot_Spline.derivative()


### first Velocity Verlet update
result_array = Velocity_Verlet(hr_vs_t[0], hv_vs_t[0], mu, RHF_Harm_Force, dt)
### first analytic result
ar_vs_t[0] = harmonic_position(om, Amp, np.pi/4, RHF_Req, 0)
### do the update N_update-1 more times
for i in range(1,N_updates):
    ### store current time
    t_array[i] = dt*i
    ### Compute VV update
    result_array = Velocity_Verlet(result_array[0], result_array[1], mu, RHF_Harm_Force, dt)
    ### store results from VV update
    hr_vs_t[i] = result_array[0]
    hv_vs_t[i] = result_array[1]
    ### compute and store results from analytic solution
    ar_vs_t[i] = harmonic_position(om, Amp, np.pi/4, RHF_Req, dt*i)

### Plot result and compare!
plt.plot(t_array, hr_vs_t, 'red', label="Velocity Verlet")
plt.plot(t_array, ar_vs_t, 'b--', label="Analytic")
plt.legend()
plt.show()
    
```
{: .language-python}


