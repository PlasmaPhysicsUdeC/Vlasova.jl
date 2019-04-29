# TODO: This is for the previus version. UPDATE

** This file is far from complete, but for now these are some basic instructions: **

To start a new simulation
-------------------------

* Choose the parameters desired in the file "parameters.jl"
* Execute "start.sh".
  To do this, and with the terminal into the program folder, type:
  
  `$ bash start.sh`

  or, simply

  `$ ./start.sh`

* In the case that you want to start the simulation from a remote terminal (with ssh, for example), you can use the tools "nohup" and "disown" so that closing the terminal won't kill the simulation. In this case, a posible syntax would be

   `$ nohup bash start.sh & disown`

   * In all both, a file called "progress" will be placed into the data folder, with an indication of the percentage of the simulation accomplished.

Some considerations
--------------------

* The data generated will be saved into the folder "data/", under the name of the simulation chosen in "parameters.jl"

* Inside the folder of the data, the program will leave a file called "plot.jl".
  Executing this file in julia,
  
  `$ julia plot.jl`
  
  will generate a few plots from the simulation data. Although the plots will only show basic/generic information, it serves well as a starting point

* Considering the posibility of errors such as electric failure, this program backs up the simulation every 10% by default.

To restart an incomplete simulation
-----------------------------------

* With the terminal into the main program folder, type:

  `$ ./continueFromBackUp.sh "simulation name"`

  and the program will do the rest.

About the algorithm
-------------------
* This program integrates the distribution function numerically separating the vlasov equation into a set of advections, as proposed by Cheng & Knorr (1976).
  For now the simulation is using the symplectic integrator BAB (also known as Velocity Verlet) from the classifiation given by Omelyan (2003). It is pending to add support for the integrator BABAB.

* Interpolations are performed using spectral methods following a backwards semi-Lagrangian scheme.


Features
--------

* Support for an arbitrary number of species is given. In "parameters.jl" there is a vector of species in which you can add a number of species to the plasma, selecting each one's name, charge, temperature and initial distribution function. All variables must be normalized to electron's magnitudes.

* The program makes use of an anisotropic filter for unphysically high frequencies into the velocity space. This allows to overcome recurrence and to avoid numerical problems coming from the Gibbs phenomenon [See https://en.wikipedia.org/wiki/Gibbs_phenomenon]

* Parallelization may be handled into "parameters.jl", changing the variables: "OMP_THREAD_NUM" and "FFTW_THREAD_NUM", which control the number of threads to use by fourier transforms (FFTW) and space/velocity advections (Programmed in FORTRAN with OpenMP).

If multiple instances of the program are run, it is a better idea to run the simulation serially (just one thread). However, under normal circumstances, a good number of threads to use is 4 for both OMP as well as FFTW.

Dependencies
------------

This program is mainly written in **Julia using the packages LinearAlgebra, Dates, FFTW and PyPlot** ( syntax for version >= 1.0 ).

Also, some parts are written in FORTRAN 90; this parts are compiled using **gfortran** and parallelized with **OpenMP**.

I hope you enjoy it!