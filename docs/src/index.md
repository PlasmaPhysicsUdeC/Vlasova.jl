# Vlasova

> A package to solve the Vlasov-Poisson system of equations numerically, under development by Jorge Gidi for the Group of Plasma Physics of the University of Concepción, Chile.

## Introduction

Vlasova is a package for julia, written in julia, which aims to provide a set of tools to perform and examine one or two-dimensional Vlasov-Poisson simulations easily.

## Features

The Vlasova package provides tools to easily:
* Define and evolve plasmas, which can be
  * 1 or 2-dimensional
  * Multi-species

* Examine the results of a simulation generically.
* Use arbitrary symplectic **non-gradient** integrators. (gradient integrators are being implemented!)
* Perform multi-threaded simulations.
  By default Vlasova will use all the threads available to the julia process.

## Asumptions

* The plasmas are electrostatic and unmagnetized.

  Vlasova considers the Vlasov equation coupled to the Poisson equation. This means that the plasmas are always assumed to behave electrostatically, since the electric field is obtained self-consistently from the electrostatic potential held by the plasma distribution function.
  Also, no magnetic field is assumed.

* The phase space is rectangular and periodic in positions and velocities

  Vlasova utilizes Fourier Transforms in the space as well as velocity dimensions, which require periodicity.
