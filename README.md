
# enhanced Green's Function Reaction Dynamics

## a modern C++ implementation


* Copyright (C) 2015-2017 AMOLF
* Copyright (C) 2005-2008 The Molecular Sciences Institute
* Copyright (C) 2009-2010 FOM Institute AMOLF
* Copyright (C) 2008-2010 RIKEN


### About this package

This package is a modern C++ implementation of the eGFRD algorithm. 

It has some important improvements over its predecessor:

* Minimal dependencies, allowing quick getting started.
* It is very fast (3 to 4 orders of magnitude faster than previous version).
* Rewritten mostly from scratch, in modern C++, using new and efficient language features.
* Improved code-style, bug- and design fixes, easier to understand and maintain.
* Multi platform support

The previous version of the eGFRD simulator had grown in the years, with many 
different developers, into a complicated and bloated codebase. For new students 
it would take several months to get familiar with the source code, and even longer 
to get the understanding and trust to change and/or add features. Furthermore the 
project had many dependencies, where some had become obsolete and others diverged, 
leading to version compatibility problems.

It became clear that this situation was no longer maintainable on the long term. 
In 2015 an effort was started to rewrite the eGFRD simulator from scratch. Resulting 
in this package.

It utilizes the new C++ language standard (modern ISO C++11 and C++14) with its 
optimized features like move-semantics and smart-pointers. Where needed types use 
fast STL-containers leading to no dependencies on external libraries (except for 
GNU Scientific Library). Most of the Template Meta Programming (TMP) model which 
was excessively present in the previous code was deserted, leaving clean 
understandable types and faster code compilation (rebuild time went down from 
15 min to less than 2 min). 

The best C++ components were cherry picked from the old codebase, like MatrixSpace 
and EventScheduler. Remaining simulation code in Python was converted to C++. 
Software engineering principles were implemented, like splitting code into a more 
modular system. The GreenFunctions, the Simulator and Logger now have their own 
separate libraries. The Bessel look-up tables are stored in data-files and not 
in-lined in the executable. This ongoing revision process now resulting in the first 
publicly available re-release of a novel and very fast simulator. 



### Authors

Modern C++ version writen by:

* Luc Blom
* Marco Seynen
* Pieter Rein ten Wolde


based on previous work done by:

* Nils Becker
* Laurens Bossen
* Kazunari Kaizu
* Moriyoshi Koizumi
* Thomas Miedema
* Sorin Tanase-Nicola
* Thomas Sokolowski
* Koichi Takahashi
* Pieter Rein ten Wolde



### License

This package is distributed under the terms of GNU General Public License
version 2.  See LICENSE.



### Building this package

See INSTALL.



### History

Koichi Takahashi initially stated development of the code in 2005 to
implement his prototype of Greens Function Reaction Dynamics
simulation method invented by Jeroen van Zon and Pieter Rein ten Wolde
in AMOLF, Amsterdam[4].  He gave a brief invited talk about
performance evaluation and applicability of the method to yeast
pheromon response pathway (the Alpha pathway) using the prototype in
the Third Annual Alpha Project Research Symposium (June 16-27, 2005, at UC
Berkeley Art Museum).

Later, in December 2006, ten Wolde, Sorin Tanase-Nicola, and Takahashi
decided to introduce the concept called first-passage processes
inspired by a paper by Opplestrup et al.[3] to Greens Function
Reaction Dynamics to further boost the performance and accuracy of the
method.  The new method was called eGFRD (enhanced Greens Function
Reaction Dynamics).  Takahashi implemented the single-body Greens
function within the year.  Tanase-Nicola derived the two-body Greens
function, and Takahashi devised and implemented a neat yet complicated
way to efficiently evaluate the function mostly in the first half of
2007.  Takahashi also implemented the main part of the code,
asynchronous discrete-event-driven kinetic monte-carlo, which has been
mostly finished within 2007, and the dynamic switching between Greens
function and brute-force Brownian Dynamics as a means to recover from
particle squeezing conditions by April 2008.  At this point,
development and debugging of the initial version (version 0.1) of the
code had been largely finished and had been used in simulation experiments
in study of various biochemical systems including dual-phosphorylation
pathways[4].

Thomas Miedema and Laurens Bossen, while masters students in the group of
Pieter Rein ten Wolde at AMOLF, added support for reaction-diffusion on
and with 1D and 2D surfaces. Laurens implemented the 1D and 2D Green's
functions in C++, Thomas implemented the algorithm in Python.

In 2009 Thomas Sokolowski and Nils Becker joined the project. Thomas S.
will extend the scheme to be able to simulate active transport processes
via molecular motors. This requires the calculation of new Green's
functions starting from the diffusion-drift equation. Nils B. recently
started working on the interplay of DNA sliding and 3D diffusion.



### Plans

Some features planned to be added are; migration of the surfaces / interactions
code (not available yet), model-entry frontend, simulations of polymers with 
multiple binding sites and states (fold reaction-rules to overcome combinatorial 
explosion of possible speciesType).

In addition to new features we are also working on improving the codebase; 
clean-up, commenting and documenting, unit-testing. Keeping up with compiler and 
library improvements. Adding more examples and usage information.



### References

* [1] Green's-function reaction dynamics: a particle-based approach for
    simulating biochemical networks in time and space; 
    van Zon and ten Wolde, J. Chem. Phys. 123 (2005).
* [2] Simulating biochemical networks at the particle level and in
    time and space: Green's function reaction dynamics; 
    van Zon and ten Wolde, Phys. Rev. Lett. 94 (2005).
* [3] First-Passage Monte Carlo Algorithm: Diffusion without All the Hops,
    Opplestrup et al., Phys. Rev. Lett. 97 (2006).
* [4] Spatio-temporal correlations can drastically change the response of a
    MAPK pathway; Takahashi, Tanase-Nicola, ten Wolde, arXiv:0907.0514v1 (2009)
