
---
title: "_ab_ _initio_ potential energy surface"
teaching: 0
exercises: 0
questions:
- "How can we simulate the vibrational motion of a molecule?"
objectives:
- "We can solve Newton's equations of motion for relative motion of atoms subject to intramolecular forces."
keypoints:
- "We can use quantum chemistry to obtain accurate intramolecular forces, and we can use molecular dynamics simulations 
to solve Newton's equations of motion to study the vibrational motion moving according to intramolecular force."
---
FIXME

{% include links.md %}

``` {: .language-python}
import numpy as np
import psi4
from matplotlib import pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
```

