function [clusteriness, p, dists, k] = pairsClusterTest(vectors, params)
% [pairsVal, p, dists, k] = pairsClusterTest(vectors [, params])
%
% Implement the PAIRS clustering test. PAIRS stands for Projection Angle
% Index of Response Similarity.
%
% Tests whether a dataset in n dimensions is likely to be spherically
% symmetric (i.e., randomly distributed on the surface of an n-1 sphere).
% There is no really good known way to do this that I could find. Here, we
% use k-nearest neighbor angles. That is, for each datapoint, we find the
% angle to its k nearest neighbors. This yields as many angles as there
% were datapoints times k. We then generate spherically random data of the
% same size and dimensionality, and perform the same operations. This gives
% us a null distribution of k-nearest-neighbor angles. For each datapoint,
% we average its k angles, then take the median of the distribution. We can
% then simply use the random distribution to determine whether the
% distribution of angles is different for the real data than the random
% data. We compare the medians, because clustering produces a skewed
% distribution. The p-value returned is two-sided.
%
% k is selected automatically unless overridden. To do the automatic
% selection, larger and larger values of k are tried until the median angle
% is >= the parameter targetAngle (see below). Default pi/4, but smaller
% values should probably be used when working in low-D (e.g., pi/8 for
% dimensionality 2 or 3).
%
% To illustrate use, imagine we are trying to find whether the responses of
% neurons tend to fall into categories. We would first build a matrix, A,
% of size CT x N, C the number of conditions, T the number of time points,
% and N the number of neurons. We would then perform PCA down to nDims:
%
% nDims = 8;                 % for example
% coeff = princomp(A);
% coeff = coeff(:, 1:nDims);
%
% Finally, we would run PAIRS on this matrix:
%
% [pairsVal, p, dists, k] = pairsClusterTest(coeff);
% 
%
% INPUTS:
%
% vectors -- the data. Should be nDatapoints x nDims
% params  -- an optional struct. It may have fields:
%   .k           -- PAIRS uses the k nearest neighbors. If this value is
%                   specified and non-NaN, it will override the automatic
%                   selection unless .targetAngle is also specified.
%   .targetAngle -- for selection of k. Overrides direct specification of
%                   k. Will ensure that the mean angle for random vectors
%                   is >= this value. Default pi/4. When working in very
%                   low-D (e.g., 2 or 3), a smaller value is recommended
%                   (e.g., pi/8).
%   .nDims       -- if present, will truncate the data to params.nDims x
%                   nDatapoints by taking the first nDims rows
%   .nBootstraps -- number of bootstraps (really simulations) to run when
%                   generating the null distribution. Default 10000.
%   .reseedRNG   -- whether to re-seed rng and set method to 'twister'.
%                   Default 1.
%   .minLength   -- minimum magnitude of vectors to include. Default 0, but
%                   will always remove zero-length vectors. This can be
%                   used to remove datapoints that are captured poorly by
%                   the PCA space.
%   .showVectors -- whether to show a 2-D projection of the normalized
%                   vectors (first 2 dims). Default 1.
%   .showHist    -- whether to show the histogram comparing the data and
%                   the Monte Carlo simulations. Default 1.
%   .histBins    -- number of bins to use in the histogram (if applicable).
%                   Default 40.
%
% OUTPUTS:
%
% pairsVal     -- 1 if perfectly clustered, 0 if perfectly random. Note
%                 that 'perfectly clustered' means only that each point has
%                 an identically angled neighbor, not that there is only a
%                 single cluster. May be negative if data is more uniform
%                 than chance.
% p            -- p-value (two-sided).
% dists        -- a struct, contains the distributions of the
%                 nearest-neighbor angles.
%   .anglesData      -- k x nDatapoints nearest-neighbor angles for the data
%   .neighbors       -- k x nDatapoints vector of which vectors made the
%                       smallest angles with this vector.
%   .anglesBootstrap -- nBoostraps x nDatapoints, nearest-neighbor angles
%                       for each Monte Carlo run (averaged over k)
% k            -- the value of k that was used (that is, this tells you
%                 that the k nearest neighbor angles were measured).
%
% Version release date: 11/11/14
%
% Copyright (c) Matt Kaufman 2013, 2014
% Cold Spring Harbor Laboratory
% antimatt AT gmail DOT com
%
% If used in published work, please cite:
% Raposo D*, Kaufman MT*, Churchland AK (2014). "A category-free neural
% population supports evolving demands during decision-making." Nature
% Neuroscience.


% Wish to match Matlab's convention for inputs (using the coeff output of
% princomp), but the transpose is easier and faster to compute on
vectors = vectors';

if size(vectors, 2) < 2
  error('pairsClusterTest:tooFewVecs', 'Must supply at least two vectors');
end


%% Defaults

k = NaN;

targetAngle = pi / 4;

nSimsForKSearch = 10;

% Dimensionality
nDims = size(vectors, 1);

% # of bootstraps (really random simulations) to run
nBootstraps = 500;

reseedRNG = 0;

minLength = 0;

showVectors = 0;
showHist = 0;

histBins = 40;

% Assign args from params struct if present
% Note that this call will potentially overwrite the values of any of the
% parameters above
if exist('params', 'var')
  warnBadAssigns(assignFromParamsStruct(params, who));
end

% Make sure targetAngle is ok
if isnan(targetAngle) || isinf(targetAngle)
  error('pairsClusterTest:badTargetAngle', ...
    'If targetAngle is specified, it should be non-NaN and non-Inf');
end


% Trim vectors if required
if nDims < size(vectors, 1)
  vectors = vectors(1:nDims, :);
elseif nDims > size(vectors, 1)
  warning('pairsClusterTest:tooFewDims', ...
    'Dimensionality of data is lower than requested for testing, using dimensionality of data');
end


%% Normalize, screen vectors

% Get magnitudes
lengths = sqrt(sum(vectors .^ 2));

% Normalize
vectors = bsxfun(@rdivide, vectors, lengths);

% Remove vectors that had zero length or were shorter than minLength
vectors = vectors(:, ~isnan(vectors(1, :)) & lengths >= minLength);



%% Pre-allocate, initialize

nVecs = size(vectors, 2);
minAnglesBootstrap = NaN(nBootstraps, nVecs);

if reseedRNG
  rng(0, 'twister');
end



%% Display vectors

if showVectors
  blankFigure;
  hold on;
  
  xs = [zeros(1, nVecs); vectors(1, :); NaN(1, nVecs)];
  ys = [zeros(1, nVecs); vectors(2, :); NaN(1, nVecs)];
  plot(xs(:), ys(:), 'b-', 'LineWidth', 1);
  rectangle('Position', [-1 -1 2 2], 'Curvature', 1, 'EdgeColor', 'k', 'LineWidth', 1.5);
end


%% If needed, search for the correct value of k

if isnan(k)
  dotProdsAll = zeros(nVecs, nVecs, nSimsForKSearch);
  
  for s = 1:nSimsForKSearch
    % Choose spherically random vectors
    rVecs = randn(size(vectors));
    
    % Normalize
    rVecs = bsxfun(@rdivide, rVecs, sqrt(sum(rVecs .^ 2)));
    
    % Get all dot products
    dotProds = rVecs' * rVecs;
    
    % Disqualify diagonal (which will be all 1's)
    dotProdsAll(:, :, s) = dotProds - 2 * speye(nVecs);
  end
  
  % Sort dot products, largest at the top (these are the smallest angles)
  dotProdsAll = sort(dotProdsAll, 1, 'descend');
  
  % We'll try larger and larger values of k (the number of nearest
  % neighbors to average).
  for testK = 1:nVecs - 1
    testAngles = acos(dotProdsAll(1:testK, :, :));
    meanAngles = mean(testAngles, 1);
    medianAngle = median(meanAngles(:));
    
    if medianAngle >= targetAngle
      k = testK;
      break;
    end
  end
end



%% Find real distribution of angles

% Get all dot products
dotProds = vectors' * vectors;

% Disqualify diagonal (which will be all 1's)
dotProds = dotProds - 2 * speye(nVecs);

if k == 1
  [maxProdsData, maxProdsDataI] = max(dotProds);
else
  [maxProdsData, maxProdsDataI] = sort(dotProds, 1, 'descend');
  maxProdsData = maxProdsData(1:k, :);
  maxProdsDataI = maxProdsDataI(1:k, :);
end



%% Generate null distribution of angles

for b = 1:nBootstraps
  % Choose spherically random vectors
  rVecs = randn(size(vectors));
  
  % Normalize
  rVecs = bsxfun(@rdivide, rVecs, sqrt(sum(rVecs .^ 2)));
  
  % Get all dot products
  dotProds = rVecs' * rVecs;
  
  % Disqualify diagonal (which will be all 1's)
  dotProds = dotProds - 2 * speye(nVecs);
  
  % Get angles to nearest neighbors
  if k == 1
    minAnglesBootstrap(b, :) = acos(max(dotProds));
  else
    dotProdsBootstrap = sort(dotProds, 1, 'descend');
    minAnglesBootstrap(b, :) = mean(acos(dotProdsBootstrap(1:k, :)), 1);
  end
end


%% Turn dot products into angles

dists.anglesData = real(acos(maxProdsData));
dists.neighbors = maxProdsDataI;
dists.anglesBootstrap = minAnglesBootstrap;


%% Compute clusteriness

bootMedians = median(dists.anglesBootstrap, 2);
bootMean = mean(bootMedians);
dataMedian = median(mean(dists.anglesData, 1));

clusteriness = (bootMean - dataMedian) / bootMean;


%% Compute p-value using bootstrap

% Alternate ways of getting p-value:
% [~, p] = ttest2(mean(dists.anglesData, 1), dists.anglesBootstrap(:), 0.05, 'left'); % one-sided
% [~, p] = kstest2(mean(dists.anglesData, 1), dists.anglesBootstrap(:));

fracLess = mean(bootMedians <= dataMedian);
% Make the p-value two-sided
if fracLess < 0.5
  p = 2 * fracLess;
else
  p = 2 * (1 - fracLess);
end


%% Show histogram

if showHist
  nBootBins = 100;
  bins = linspace(0, pi/2, histBins);
  bootBins = linspace(0, pi/2, nBootBins);
  
  blankFigure;
  hold on;
  
  nData = histc(mean(dists.anglesData, 1), bins);
  nBoot = histc(dists.anglesBootstrap(:), bootBins) / nBootstraps * nBootBins / histBins;
  hb = bar(bins + diff(bins(1:2))/2, nData, 'BarWidth', 0.8, 'FaceColor', 'k', 'LineStyle', 'none');
  set(get(hb, 'BaseLine'), 'LineStyle', 'none');
  plot(bootBins(1:end-1) + diff(bins(1:2))/2, nBoot(1:end-1), 'r-', 'LineWidth', 2);
  
  theMax = max([nData(:); nBoot(:)]);
  theMax = ceil(theMax);
  
  axParams.fontSize = 11;
  axParams.tickLocations = [0 pi/4 pi/2];
  axParams.tickLabelLocations = axParams.tickLocations;
  axParams.tickLabels = {'0', 'pi/4', 'pi/2'};
  axParams.axisLabel = 'NN angle';
  AxisMMC(0, pi/2, axParams);
  
  axParams.axisOrientation = 'v';
  axParams.tickLocations = [0 theMax];
  axParams.tickLabelLocations = axParams.tickLocations;
  axParams.tickLabels = {'0', num2str(theMax)};
  axParams.axisLabel = 'Count';
  axParams.axisOffset = -pi/40;
  AxisMMC(0, theMax, axParams);
  
  axis tight
end