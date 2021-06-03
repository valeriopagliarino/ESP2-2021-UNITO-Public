# -*- coding: utf-8 -*-
"""Error propagation interferometer ESP2.ipynb
    
Original file is located at
    https://colab.research.google.com/drive/1nWGtJxMM0yOAu4G5nh81BKLvlYreeSzu
    
"""

import sympy as sy
import numpy as np
pi = np.pi

from sympy import *
angle_mean, nFrGlass, refLambda, ppTichness = symbols("angle_mean nFrGlass refLambda ppTichness")
angle_mean_err, nFrGlass_err, refLambda_err, ppTichness_err = symbols("angle_mean_err nFrGlass_err refLambda_err ppTichness_err")

formula = ((ppTichness-refLambda*nFrGlass/2.)*(1-cos(angle_mean*pi/180.)))/(ppTichness*(1-cos(angle_mean*pi/180.))-(refLambda*nFrGlass/2.));

e1 = (diff(formula, angle_mean) ** 2) * (angle_mean_err ** 2)
e1

e2 = (diff(formula, nFrGlass) ** 2) * (nFrGlass_err ** 2)
e2

e3 = (diff(formula, refLambda) ** 2) * (refLambda_err ** 2)
e3

e4 = (diff(formula, ppTichness) ** 2) * (ppTichness_err ** 2)
e4

err = sqrt(e1 + e2 + e3 + e4)

err

err.evalf(subs={ppTichness: 5.38e-3, ppTichness_err: 0.05e-3, refLambda: 633e-9, refLambda_err: 1e-9, nFrGlass: 50, nFrGlass_err: 0.5, angle_mean: 7.24, angle_mean_err: 0.1})

formula.evalf(subs={ppTichness: 5.38e-3, refLambda: 633e-9, nFrGlass: 50, angle_mean: 7.24})

