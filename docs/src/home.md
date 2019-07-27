# Vlasova

> A package to solve the Vlasov-Poisson system of equations numerically, under development by Jorge Gidi for the Group of Plasma Physics of the University of Concepción, Chile.

## Introduction

Vlasova is a package *mostly* written in julia, which aims to provide a set of tools to perform and examine Vlasov-Poisson simulations easily.

## Features

The Vlasova package provides tools to easily:
* Define and evolve plasmas, which can be
  * 1 or 2-dimensional
  * Multi-specie

* Examine the results of a simulation generically
* Implement *non-gradient* integrators
* Perform multi-threaded simulations

## Asumptions

* The plasmas are electrostatic.

Vlasova considers the Vlasov equation coupled to the Poisson equation. This means that the plasmas are always assumed to behave electrostatically, since the electric field is obtained self-consistently from the electrostatic potential held by plasma *distribution function*.

* The phase space is periodic in positions and velocities

Vlasova utilizes Fourier Transforms in the space as well as velocity dimensions, which require periodicity along them.