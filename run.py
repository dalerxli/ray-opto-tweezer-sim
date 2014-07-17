#!/usr/bin/env python

import subprocess
from math import pi

exec_path = ["./app"]

# Here come the arguments for the simulation

# The lens focal distance (in m)
df = 30e-3;

# Incoming beam width
l_r = 1.1e-3;

# Medium refractive index
ne = 1;

# Relative index of the sphere
sph_n = 1.5;

# Quantity of steps to divide the r and theta (integration variables) into.
r_steps = 300;
th_steps = 300;

# The limits of integration for x (radial position of particle) IN TERMS OF
# PARTICLE RADIUS. The particle position is relative to the focus.
x_init = 0;
x_final = 0;
x_steps = 100;

# The limits of integration for y (radial position of particle) IN TERMS OF
# PARTICLE RADIUS. The particle position is relative to the focus.
y_init = 0;
y_final = 0;
y_steps = 3;

# The limits of integration for z (axial position of particle) IN TERMS OF
# PARTICLE RADIUS. The particle position is relative to the focus.
z_init = -1000e-6;
z_final = 1000e-6;
z_steps = 1003;

# Double trap?
double_trap = 0

# Sphere radiuses
s_rl = [4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 10e-6, 11e-6, 12e-6]

# Light wavelength
lam = 532e-9 

# List of arguments

cmd_l = []

for s_r in s_rl:
    args_list = [
    df,
    l_r,
    ne,

    sph_n,

    r_steps,
    th_steps,

    x_init,
    x_final,
    x_steps,

    y_init,
    y_final,
    y_steps,

    z_init,
    z_final,
    z_steps,

    double_trap,

    s_r,
    lam,
    ]

    args_list = map(lambda x: str(x), args_list)
    cmd = "echo -e \'\\n\\n\"Radius {0:.1f} um\"\' >> mul2.tsv && ".format(s_r*1e6)
    cmd += " ".join(exec_path+args_list)
    cmd += " >> mul2.tsv"
    cmd_l.append(cmd)
print " && ".join(cmd_l)
