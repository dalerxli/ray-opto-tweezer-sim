#!/usr/bin/env python

import subprocess as sp
from math import pi

# Path to the simulation executable
exec_path = "./app"

# Output file (to write the simulation results)
out_path = "output.tsv"

# Here come the arguments for the simulation

# The lens focal distance (in m)
df = 30e-3;

# Incoming beam radius
l_r = 1.1e-3;

# External medium refractive index
ne = 1;

# Relative index of the sphere
sph_n = 1.5;

# Quantity of steps to divide the r and theta (integration variables) into.
r_steps = 300;
th_steps = 300;

# The limits of integration for x (radial position of particle) in meters. 
# The particle position is relative to the focus.
x_init = 0;
x_final = 0;
x_steps = 100;

# The limits of integration for y (radial position of particle) in meters. 
# The particle position is relative to the focus.
y_init = 0;
y_final = 0;
y_steps = 3;

# The limits of integration for z (axial position of particle) in meters. 
# The particle position is relative to the focus.
z_init = -1000e-6;
z_final = 1000e-6;
z_steps = 1003;

# Double (contrapropagating) trap? 1 for yes, 0 for no
double_trap = 0

# Sphere radiuses to be simulated (in meters)
s_rl = [4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 10e-6, 11e-6, 12e-6]

# Light wavelength in meters
lam = 532e-9 

# Generate the commands
for s_r in s_rl:
    # Create the argument list for subprocess, including the executable itself
    args_list = [
    exec_path,
    
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
    
    out_file = open(out_path, "a")

    args_list = map(lambda x: str(x), args_list)
    
    # Write the current particle radius for plotting software
    out_file.write('"Radius {0:.1f} um"\n'.format(s_r*1e6))
    
    # Make sure that the above is written
    out_file.flush()
    
    # Call the executable and wait for it to complete
    sp.call(args_list, stdout=out_file)
    
    out_file.flush()
    out_file.write("\n\n")
    
    print "Particle radius {0:.1f} um finished".format(s_r*1e6)
    
    # Close the file to write changes
    out_file.close()
