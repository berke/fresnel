#!/usr/bin/env python3

import numpy as np
import scipy.special

fd = open("src/test_data.rs","w")
fd.write("pub const TEST_DATA : &[[f64;3]] = &[");

def rndfloat():
    e = np.random.randint(-20,20)
    s = 2*np.random.randint(0,2)-1
    f = np.random.random_sample()
    x = s*f*(2**e)
    return x

xs = [0.0,-1.0,1.0,0.5,-0.5,2.0,-2.0]
xs.extend([rndfloat() for i in range(10000)])
print(f"xs = {xs}");

for x in xs:
    (c,s) = scipy.special.fresnel(x)
    fd.write(f"[{x:.20},{c:.20},{s:.20}],\n")

fd.write("];");
