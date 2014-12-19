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
The various tests can be accessed through three interfaces, `DepTest1`, 
`DepTest2` and `UniSphereTest` for one-sample tests, two-sample tests and 
one-sample tests one the sphere, respectively.

## Multivariate (In)dependence, Sphericity and Homogeneity
```
% Independent, but non-spherical data
sigma = diag([ones(1,25),0.5*ones(1,5)]);
x = (sigma*randn(50,30)')';

% Independence tests
DepTest1(x,'test','spearman')
DepTest1(x,'test','kendall')

% Sphericity tests
DepTest1(x,'test','john')
DepTest1(x,'test','wang')
DepTest1(x,'test','sign')
DepTest1(x,'test','bcs')
```

```
% Non-indepedent data, with ~0 correlation, same distribution
x = rand(200,1); y = rand(200,1);
xx = 0.5*(x+y)-0.5; yy = 0.5*(x-y);
corr(xx,yy)

% Two-sample Independence tests
DepTest2(xx,yy,'test','dcorr')
DepTest2(xx,yy,'test','hsic')

% Do the samples come from the same distribution?
DepTest2(xx,yy,'test','mmd')
DepTest2(xx,yy,'test','energy')
```

```
% Independent data, different distributions
x = randn(200,1); y = rand(200,1);

% Two-sample Independence tests
DepTest2(x,y,'test','dcorr')
DepTest2(x,y,'test','hsic')

% Do the samples come from the same distribution?
DepTest2(x,y,'test','mmd')
DepTest2(x,y,'test','energy')
```

## Uniformity on hypersphere
```
% Non-uniform samples
sigma = diag([1 5 1]);
x = (sigma*randn(50,3)')';

UniSphereTest(x,'rayleigh') % Note failure of Rayleigh test, since resultant is zero
UniSphereTest(x,'gine-ajne') 
UniSphereTest(x,'randproj') 
UniSphereTest(x,'bingham') 
```

Contributions
--------------------------------
Copyright (c) 2014 Brian Lau [brian.lau@upmc.fr](mailto:brian.lau@upmc.fr), see [LICENSE](https://github.com/brian-lau/highdim/blob/master/LICENSE)

Please feel free to [fork](https://github.com/brian-lau/highdim/fork) and contribute!