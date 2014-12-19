highdim
==========

A Matlab library for statistical testing with high-dimensional data, including
one and two-sample tests for homogeneity, uniformity, sphericity and
independence. Of note are implementations of a number of modern methods 
appropriate for dealing with data where dimensionality is large, possibly
exceeding the number of samples.

# Installation
Download [highdim](https://github.com/brian-lau/highdim/archive/master.zip), 
add the resulting folder to your Matlab path, and you're ready to go. 
Folders with a `+` are packages that do not need to be explicitly added to your path, 
although its [parent folder does](http://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html#brfynt_-3).

# Examples
The various tests can be accessed through three interfaces, `DepTest1`, `DepTest2` and `UniSphereTest` for one-sample tests, two-sample tests and tests one-sample tests one the sphere, respectively.
## Multivariate (In)dependence

## Multivariate Homogeneity

## Sphericity

## Uniformity on hypersphere


Contributions
--------------------------------
Copyright (c) 2014 Brian Lau [brian.lau@upmc.fr](mailto:brian.lau@upmc.fr), see [LICENSE](https://github.com/brian-lau/highdim/blob/master/LICENSE)

Please feel free to [fork](https://github.com/brian-lau/highdim/fork) and contribute!