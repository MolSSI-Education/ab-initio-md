---
title: "Running the molecular dynamics simulation!"
teaching: 0
exercises: 0
questions:
- "What does the molecular dynamics simulation tell us?"
objectives:
- "To demonstrate the execution of a molecular dynamics simulation and analysis of basic results."
keypoints:
- "The direct output of a molecular dynamics simulation are the positions and velocities/momenta of the atoms being studied."
---

Finally, let's simulate the vibrational motion of HF subject to the *ab* *initio* forces we computed earlier and compare them to the Harmonic motion; recall we have already obtained spline objects for RHF, MP2, and CCSD forces called RHF_Force, MP2_Force, and CCSD_Force!
We will also initialize the simulations using the same values as we did with the Harmonic case to aid our comparison.

```
### how many updates do you want to perform?
N_updates = 200000

### establish time-step for integration to be 0.02 atomic units... this is about 0.0005 femtoseconds
### so total time is 200000*0.02 atomic units of time which is ~9.6e-13 s, or 960 fs
dt = 0.02

### Now use r_init and v_init and run velocity verlet update N_updates times, plot results
### these arrays will store the time, the position vs time, and the velocity vs time
r_vs_t = np.zeros(N_updates)
v_vs_t = np.zeros(N_updates)
t_array = np.zeros(N_updates)

### harmonic results
hr_vs_t = np.zeros(N_updates)
hv_vs_t = np.zeros(N_updates)


### first entry is the intial position and velocity
r_vs_t[0] = r_init
v_vs_t[0] = v_init

hr_vs_t[0] = r_init
hv_vs_t[0] = v_init

### first Velocity Verlet update
result_array = Velocity_Verlet(r_init, v_init, mu, RHF_Force, dt)
#hresult_array = Velocity_Verlet(r_init, v_init, mu, HF, dt)

print(result_array)
#print(hresult_array)

### do the update N_update-1 more times
for i in range(1,N_updates):
    tmp = Velocity_Verlet(result_array[0], result_array[1], mu, RHF_Force, dt)
    result_array = tmp
    t_array[i] = dt*i
    r_vs_t[i] = result_array[0]
    v_vs_t[i] = result_array[1]
    #tmp = Velocity_Verlet(hresult_array[0], hresult_array[1], mu, HF, dt)
    #hr_vs_t[i] = hresult_array[0]
    #hv_vs_t[i] = hresult_array[1]

### Plot the trajectory of bondlength vs time:
#plt.plot(t_array, r_vs_t, 'red', t_array, hr_vs_t, 'blue')
plt.plot(t_array, r_vs_t)
plt.show()
```
{: .language-python}
