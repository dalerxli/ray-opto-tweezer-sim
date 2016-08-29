# ray-opto-tweezer-sim
##Description
`ray-opto-tweezer-sim` is a program written for simulating forces on homogeneous dielectric spheres trapped in optical tweezers where the ray optics approximation is applicable, i.e. when the particle is significantly bigger than the wavelength of the trapping light and is bigger that the waist of the light beam.

This code is written in Python3+NumPy, which should make it easy to install and use (see "Installation" below), while also being reasonably fast (courtesy of NumPy): evaluating the force on a particle is done at 20 particle positions per second on a single core of Intel Celeron 1007U (an entry-level laptop CPU).

The theoretical background for this code is described in Ashkin, 1992. Since I couldn't find the code used by Ashkin in his paper, I have decided to write this one.

### Current features
- Calculation of the Q-factor (see Ashkin, 1992 for the definition) is done at a relatively high speed and with controllable precision. The Q factor is computed relatively to the power that passes through the lens (not the total incident power on the lens), as this approach seems to be more convenient.

- The program (now) accepts a TEM00 Gaussian beam with an arbitrary beam waist specified by the ratio between the beam waist and the lens radius. That way, by specifying the NA of the focusing lens, the cone of light is completely determined. In case if the Gaussian beam is significantly larger than the lens, it is clipped. That way, a uniform lens illumination can be considered by setting a to a large value (say, 1000), so that uniform illumination cases are also considered.

- Any spatially uniform polarization state of the incoming light is possible, such as linear polarization (e.g. *(1,0)* in Jones notation), circular (equivalent here to non-polarized), specified by *(1,1j)* where *j* is the imaginary unit. More general elliptic polarization states are also possible (equivalent here to arbitrary mixtures of p and s polarizations).

- Easy setting-up: all the simulation values and the output file are set in the amply-commented configuration file "config.ini".

- Particle positions to be evaluated are specified as ranges in each of the X,Y,Z coordinates, which allows to construct 1D, 2D and 3D force profiles with an arbitrary number of steps in each of the coordinates, while increasing the performance (some calculations can be recycled by the program).

- Trustworthy calculations: an automated test suite (in `tests` subdirectory) verifies the operation of the code, beginning with the basics (Snell law and Fresnel coefficients implementation) and ending by the total force exerted on a particle. Where possible, the test suite ensures that the generated results are consistent with the published ones (again, with Ashkin, 1992). Special attention is paid to corner cases for each of the simulation functions in order to guarantee their correctness.

### Upcoming features
- Arbitrary beam intensity profile that can be easily specified by the user.

- Arbitrary spatial polarization profile to allow simulating e.g. radial polarization.

- Multithreading support to speed things up on multi-core or multi-CPU machines.

- Arbitrary ray-transfer function for the lens to allow simulating e.g. lenses with aberrations or lenses more interesting than the occasional Abbe-compliant objectives.

## Installation
### Installation of requirements
#### Universal way

In order to use this program, Python 3, together with some additional packages (Scipy and NumExpr), is required. The easiest way to install it (Windows, Linux or Mac OS X) is to install the [Anaconda](https://www.continuum.io/downloads) Python distribution that already includes all the necessary packages. Just make sure to download the Python 3 version of Anaconda as opposed to Python 2.

#### Linux
Alternatively, users of Linux can easily install Python 3 and the necessary packages by using their distribution's package system. In Ubuntu (at least in 14.04), for example, all the necessary packages can be easily installed with the console command:

```
sudo apt-get install python3-numexpr python3-scipy
```

This command will also install all the dependencies and, if it is not installed yet, Python3 itself. This approach has the advantage of a tight integration with the operating system and may not require to download as much data.

### Setting up of the program
The program requires almost no preparations to be run for the first time. The only thing you have to do is to download the code itself by pressing the "Clone or download" green button on this page, and then selecting "Download ZIP". After that, just extract the contents of the ZIP file anywhere you'd like to.

If you have used Git, then, equivalently, you can also clone this repository (I assume you know how to do that). Also, cheers!

### Testing the program
The program comes configured for a first run, which is useful to check that it is working fine on your machine. In this run, it calculates 100 positions of a particle along the beam axis (Z) (from -2 to 2, where the position is measured in terms of the particle radius). This calculation only takes 5 seconds even on my relatively weak laptop, so it shouldn't take a lot of time.

To perform this first run, first navigate into the directory in which you extracted the downloaded ZIP file. If you don't know how, then simply open a console or a command prompt, and execute a command like that:

```cd directory/you/are/looking/for```

In case of Windows, replace "/" into "\". The directory you should navigate to should contain the "run.py" file.

After that, run the following command:

```python3 run.sh```

If the command ran without any errors, then everything is fine (some warnings *may* appear).

If you see "No module named ..." errors, then check that you have installed your Python modules correctly. On Anaconda, this should not happen. If other errors appear, then please raise an issue in this repository ("Issues" tab on the top of this page).

After a successful run, a "results.tsv" file will appear that contains the data. The first three columns are the coordinates of the particle and the next three are the force acting on the particle in this position. You can open it in nearly any software that can plot graphs (e.g. Excel, gnuplot, Mathematica and so on).

## Usage
Now that the program works, you can modify the parameters you would like to in the file "config.ini". This file describes every parameter with detail, so just open it and have fun. After modifying "config.ini", just run `python3 run.sh` to generate data according to the new configuration. Careful to not overwrite previous data (you can change the output file in "config.ini").