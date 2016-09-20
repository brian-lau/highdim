function [aInd, p, varToThisDimInTSpace, varToThisDimInSSpace, varToThisDimWorstCase, randCI, sigs] = ...
  varAlignment(ASpace, ATest, nDims, showStars, makePlot, nRandSeeds)
% [alignInd, p, varToThisDimTest, varToThisDimSpace, varToThisDimWorstCase, randCI, sigs] = ...
%   varAlignment(ASpace, ATest, nDims [,showStars] [, makePlot] [, nRandSeeds])
%
% This function is looking to see whether the data in ATest lives in the
% same space as the data in ASpace. It does so by implementing Variance
% Alignment analysis from Raposo*, Kaufman*, and Churchland 2014.
%
% ASpace and ATest should be of size nObservations x nVariables. nDims
% determines the dimensionality of the initial PCA step (see below).
%
% This is a tricky thing to do. If, for example, one simply ran PCA on both
% and compared the angles between the first PC of each space, in a high-D
% space that would tend to be 90 degrees (especially if a few dimensions
% have nearly equal variance). Moreover, what space the data lives in is
% not a clean thing; it's a question of how much variance is in each
% dimension.
%
% Here is the key logic of Variance Alignment. If the two spaces are
% aligned, we expect that the dimensions that capture the most variance for
% one space will capture a lot of variance for the other. If they are
% actively misaligned, we expect that the best dimensions for one space
% will be the worst dimensions for the other.
%
% To follow through on this intuition, this function does the following:
%
% 1) PCA the A matrix containing both 'test' and 'space' data, down to
%    nDims dimensions.
% 2) In this PCA'd space, PCA just the 'space' data to get the orientation
%    of its variance ellipse.
% 3) Rotate a copy of the 'test' data into this 'space' space.
% 4) Rotate a copy of the 'test' data into the test-only PCA space.
% 5) for d = 1:nDims, check how much variance is in the test data in the
%    'space' space up to d dimensions divided by how much variance is in
%    'test' space up to d dimensions.
%
% This fraction is plotted as a function of d, as a thick black line in the
% plot.
%
% For perfectly aligned spaces, the fraction above should always be 1
% (black dashed line in the plot). For maximally misaligned spaces, this
% should be the same as if in step 4 we used the 'test' space and ordered
% the dimensions from lowest var to highest instead of highest to lowest
% (blue line in the plot). For comparison, we can also choose random
% spaces. The mean of many random choices is shown in thick red, with 95%
% CIs (not corrected for multiple comparisons) as a gray fill. Note that
% for all space choices, the fraction will be bounded by [0, 1], and the
% value at d == nDims must always be unity, since we have run out of
% dimensions for variance to hide in.
%
% We also derive an overall index of space alignment. +1 means perfectly
% aligned, -1 means maximally misaligned, 0 means random. To compute this
% index, we subtract the random line from the data line, and average over
% the resulting values. If the average is positive, we divide by [(perfect
% alignment) - (the random line)]. If the result was negative, we divide by
% [(random) - (maximally misaligned)] . We can get the random distribution
% of the index to compute a p value, which is two-sided.
%
% To test each individual point for significance we use a two-sided test at
% an alpha level of 0.05, using a bootstrap (10,000 iterations by default)
% and the Holm-Sidak procedure. This corrects for multiple comparisons, and
% is strictly better than Bonferroni or Holm-Bonferroni. Significant points
% are shown with a star on the data line, if showStars is 1 (default 1).
% The last point is not tested, because it will always be unity, and is
% unstarred.
%
% Version release date: 11/11/14
%
% Copyright (c) Matt Kaufman 2013, 2014
% Cold Spring Harbor Laboratory
%
% If used in published work, please cite:
% Raposo D*, Kaufman MT*, Churchland AK (2014). "A category-free neural
% population supports evolving demands during decision-making." Nature
% Neuroscience.

%% Error check

if nargin < 3
  error('varAlignment:tooFewArguments', 'Too few arguments to varAlignment');
end

if size(ASpace, 2) ~= size(ATest, 2)
  error('varAlignment:variableMismatch', 'Number of columns in ASpace and ATest must match');
end

if size(ASpace, 1) < size(ASpace, 2) || size(ATest, 1) < size(ATest, 2)
  warning('varAlignment:badSize', ...
    'The A matrices are probably transposed. Should be observations x variables');
end


%% Optional arguments

if ~exist('makePlot', 'var')
  makePlot = 1;
end

if ~exist('nRandSeeds', 'var')
  nRandSeeds = 1000;
end

if ~exist('showStars', 'var')
  showStars = 1;
end



%% Parameters

alph = 0.05;  % Statistical significance alpha

gray = 0.85;  % Gray level for displaying chance fills, higher is lighter


%% Do initial PCA on combined A matrix to get a starting space, break back up

[~, scores] = princomp([ASpace; ATest]);

ASpace = scores(1:size(ASpace, 1), 1:nDims);
ATest = scores(size(ASpace, 1)+1 : end, 1:nDims);


%% Project test data into the spaces found by PCA'ing either test or 'space' data
% NOTE: all the 'var' variable below are examining the "test" data.

[~, testInTest] = princomp(ATest);

coeffSpace = princomp(ASpace);
testInSpace = ATest * coeffSpace;


%% Get denominator, of maximally aligned case
varToThisDimInTSpace = cumsum(var(testInTest));


%% What's in the data
% NOTE: this is examining the "test" data, but in the "space" space. Making
% this clearer made the variable name long.
varToThisDimInSSpace = cumsum(var(testInSpace));

dataLine = varToThisDimInSSpace ./ varToThisDimInTSpace;


%% Most misaligned case
varToThisDimWorstCase = cumsum(fliplr(var(testInTest)));


%% Random case

%rng(0, 'twister');

C = cov(ATest);

randomCase = zeros(nRandSeeds, length(varToThisDimInTSpace));
for r = 1:nRandSeeds
  % Choose random coeffs
  coeffRand = orth(rand(nDims));
  
  % Get variance along those dimensions
  randomCase(r, :) = diag(coeffRand' * C * coeffRand);
end

varToThisDimRand = cumsum(mean(randomCase));

% Could correct quantiles for multiple comparisons as below, but it's not
% clear that this makes more sense. Instead, just correct when testing for
% significance
% Other versions:
% quantCMC = alph / (2 * nDims);   % Bonferroni correction
% quantCMC = 1 - (1 - alph/2) ^ (1 / nDims);  % Sidak correction
% randCI = quantile(cumsum(randomCase, 2), [quantCMC 1-quantCMC]);
randCum = cumsum(randomCase, 2);
randCI = quantile(randCum, [alph/2 1-alph/2]);


%% Determine significance

% Note: using non-normalized values for both the data and the random cases.
% Equivalent to using normalized for both, but faster.
% Also, we'll exclude the last value, which is always 1.
sigs = holmSidakOnBootstrap(varToThisDimInSSpace(1:end-1), randCum(:, 1:end-1), alph);
sigs = [sigs false];



%% Plot, if requested

if 0%makePlot
  blankFigure;
  hold on;
  
  % Random case CIs
  fill([1:nDims nDims:-1:1], ...
    [randCI(1, :)./varToThisDimInTSpace fliplr(randCI(2, :)./varToThisDimInTSpace)], ...
    gray*[1 1 1], 'LineStyle', 'none');
  
  % Ceiling (unity)
  plot(ones(1, length(varToThisDimInTSpace)), 'k--');
  
  % Most misaligned case
  maxOrth = varToThisDimWorstCase ./ varToThisDimInTSpace;
  plot(maxOrth, 'b-');
  
  % Random case
  randLine = varToThisDimRand ./ varToThisDimInTSpace;
  plot(randLine, 'r-', 'LineWidth', 2);
  
  % What's in the data
  plot(dataLine, 'k-', 'LineWidth', 2);
  
  % Make plot nice
  addAxes(nDims);
  
  ylim([-0.03 1.1]);
  
  applyLabels(dataLine, randLine, maxOrth, randCI(2, :) ./ varToThisDimInTSpace, gray);
  
  if showStars
    drawStatStars(sigs, dataLine);
  end
else
    maxOrth = varToThisDimWorstCase ./ varToThisDimInTSpace;
  randLine = varToThisDimRand ./ varToThisDimInTSpace;

end


%% Compute index and index stats

aInd = getAlignIndex(dataLine, randLine, maxOrth);
fprintf('Alignment index: %0.3f\n', aInd);

% Index stats
randLines = bsxfun(@rdivide, randCum, varToThisDimInTSpace);
randInds = getAlignIndex(randLines, randLine, maxOrth);
p = mean(randInds > aInd);
if aInd < 0
  p = 1 - p;
end
p = p * 2;

if makePlot
  fprintf('p (two sided) = %0.5f\n', p);
end



function addAxes(nDims)

% Horizontal
axParams.tickLocations = 1:nDims;
axParams.extraLength = 1.5;
axParams.fontSize = 16;
axParams.axisLabel = 'Cumulative dimensions';
axParams.axisOffset = -0.006;

AxisMMC(1, nDims, axParams);

% Vertical
axParams.axisOrientation = 'v';
axParams.tickLocations = [0 0.5 1];
axParams.axisLabel = 'Var to this dim in train space / in test space';
axParams.axisLabelOffset = (nDims - 1) * 0.08;
axParams.axisOffset = 1 - (nDims - 1) * 0.02;

AxisMMC(0, 1, axParams);



function drawStatStars(sigs, dataLine)

for d = find(sigs)
  plot(d, dataLine(d), 'k*', 'MarkerSize', 10, 'LineWidth', 2);
end


function applyLabels(dataLine, randLine, maxOrth, randCITop, CIColor)

nDims = length(dataLine);

x = 1 + 0.95 * (nDims - 1);
text(x, 1.03, 'Perfect alignment', 'HorizontalAlignment', 'right');

x = 1 + 0.25 * (nDims - 1);
text(x, 0.03 + valAtX(x, dataLine), 'Data', 'HorizontalAlignment', 'right', 'FontSize', 16, 'FontWeight', 'bold');

x = 1 + 0.4 * (nDims - 1);
text(x, 0.03 + valAtX(x, randLine), 'Random', 'Color', 'r', 'HorizontalAlignment', 'right');

x = 1 + 0.7 * (nDims - 1);
text(x, -0.03 + valAtX(x, maxOrth), sprintf('Maximally\nmisaligned'), 'Color', 'b', 'VerticalAlignment', 'top');

x = 1 + 0.6 * (nDims - 1);
text(x, -0.05 + valAtX(x, randCITop), 'Random 95% CIs', 'Color', CIColor/3*[1 1 1], 'HorizontalAlignment', 'right');



function y = valAtX(x, ys)
% Find y value on a linearly interpolated function for a non-integer value
% of x

belowX = floor(x);
aboveX = ceil(x);

belowY = ys(belowX);
aboveY = ys(aboveX);

if belowY == aboveY
  y = belowY;
else
  y = belowY + (x - belowX) * (aboveY - belowY) / (aboveX - belowX);
end



function sigs = holmSidakOnBootstrap(dataVals, randVals, alph)
% sigs = holmSidakOnBootstrap(dataVals, randVals, alph)
%
% Perform Holm's procedure over Sidak corrections for multiple comparisons

nDims = length(dataVals);

sigs = false(1, nDims);

while ~all(sigs)
  testCols = ~sigs;
  
  n = sum(testCols);
  alphaCMC = 1 - (1 - alph/2) ^ (1 / n);
  
  randCI = quantile(randVals(:, testCols), [alphaCMC 1-alphaCMC]');
  
  % Note: to test for being outside the confidence intervals, enforce
  % non-zero margin in case confidence intervals are equal.
  newSigs = dataVals(testCols) < randCI(1, :) - 1e-14 | dataVals(testCols) > randCI(2, :) + 1e-14;
  sigs(testCols) = newSigs;
  
  if ~any(newSigs)
    return;
  end
end



function aInd = getAlignIndex(dataLines, randLine, maxOrth)

rawVal = sum(bsxfun(@minus, dataLines, randLine), 2);

aInd = zeros(1, size(dataLines, 1));

pos = rawVal > 0;
aInd(pos) = rawVal(pos) / sum(1 - randLine);

neg = rawVal < 0;
aInd(neg) = rawVal(neg) / sum(randLine - maxOrth);
