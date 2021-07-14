This repository contains some template code that may be useful for algorithms-related exercises from the July 2021 Bad Honnef Physics School, [Methods of Effective Field Theory and Lattice Field Theory](https://www.dpg-physik.de/veranstaltungen/2020/methods-of-effective-field-theory-and-lattice-field-theory).  There's no need to use it if you already have a preferred codebase.  Eventually model solutions will be added in a separate branch.

## Hybrid Monte Carlo
The goal is to write and test a hybrid Monte Carlo (HMC) algorithm for scalar phi^4 field theory in three dimensions.  We can break this up into several small exercises, adapting larger-scale [tutorials](https://nic.desy.de/sites2009/site_nic/content/e44192/e62778/e91179/e91180/hmc_tutorial_eng.pdf) run by Stefan Schaefer at the 2009 Les Houches Summer School, [Modern perspectives in lattice QCD](https://nic.desy.de/e44192/e62778).

First are fundamental coding tasks:
* Generate gaussian-distributed momenta
* Implement the computation of the effective hamiltonian `H`
* Derive the 'force' needed to update the momenta
* Implement the elementary updates of the scalar field variables and their conjugate momenta
* Combine the elementary updates into a leap-frog integrator
* Implement the computation of `ΔH` and set up the accept/reject test

Now we can test the code (optionally using the model solution to appear):
* Check the reversibility of the integrator
* Check that the root-mean-square `ΔH` scales appropriately with the step size
* Check the Creutz equality
* Combine the elementary updates into an Omelyan--Mryglod--Folk integrator and compare the root-mean-square `ΔH` vs. the leap-frog integrator

Once the tests are passed, we can compute more interesting things:
* Measure the magnetization
* Monitor the equilibration of the magnetization across the critical line
* Estimate the auto-correlation time of the magnetization across the critical line
* Experiment as curiosity drives you ([hep-lat/9806012](https://arxiv.org/abs/hep-lat/9806012) might provide inspiration)

## Conjugate gradient [pending]
The goal is to write a conjugate gradient (CG) algorithm to invert a lattice fermion operator.  Since lattice fermions are trickier than scalar field theories, a more complete framework is provided for those who don't have an existing codebase they would prefer to work with.

The provided code considers the tight-binding hamiltonian of graphene on the spatial honeycomb lattice.  It is not suitable for serious research, in part because the fermion operator is implemented as a full matrix rather than a function efficiently implementing the matrix--vector operation.  On the other hand, this may make the CG implementation more transparent.
