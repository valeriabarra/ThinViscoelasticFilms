# Code for Thin Viscoelastic Film simulations on flat substrates

This repository contains a Fortran 90 code for the numerical simulations of thin viscoelastic films and drops on solid, flat substrates. The nomenclature used in the code, reflects the reference paper 
[Barra et al., "Interfacial dynamics of thin viscoelastic films and drops", JNNFM (2016)](http://dx.doi.org/10.1016/j.jnnfm.2016.10.001). 
All terms and nomenclature in the code refer to that publication. Please ask the author if you cannot access the reference paper.

This code is based on the finite difference discretization, on a staggered grid. 
It uses a fixed spatial grid size $\Delta x$ and an adaptive time step $\Delta t$.

Different contributors during the years have worked on different parts of this code. We acnlowledge:

Te-Sheng Lin

Nebo Murisic

Ivana Seric

Michael Lam

Valeria Barra

Ryan Allaire

all supervised or co-supervised by Prof. Lou Kondic, 
from the Department of Mathematical Sciences, at the New Jersey Institute of Technology.

This code solves the thin film (or long wave) approximation of thin viscoelastic films 
on a solid, flat substrate subject to the van der Waals interaction force, in two spatial dimensions. 
The governing equations are obtained within a long-wave (lubrication) approximation of the Navier–Stokes equations 
with Jeffreys model for viscoelastic stresses. This version of the code does not include any gravitational effects. 

The implementation provided is suitable for CPUs in serial mode. 

To build the prject, move into the [Debug](./Debug/) folder by running 

`cd Debug`

and then either run the provided compile script, by running

`bash CompileScript`

or simply build with 

`make`

To execute the program, run with the name of the executable created

`./ThinVEFilm`

or launch in a scheduling system.

The output files are intended to be placed in an `Output` folder, on the same level of the `Debug` folder.


Users and contributors are always welcome.
