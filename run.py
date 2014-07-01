#!/usr/bin/env python

import subprocess

exec_path = ["./app"]

# Here come the arguments for the simulation

# The lens focal distance (in m)
df = 30e-3;

# Incoming beam width
l_r = 3e-3;

# Medium refractive index
ne = 1;

# Relative index of the sphere
sph_n = 1.33;

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

# Sphere radius
s_r = 7e-6

# Light wavelength
lam = 532e-9 

# List of arguments
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
print " ".join(exec_path+args_list)
