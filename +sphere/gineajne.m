% GINEANGE                    Sum of Gine, Anje statistics for spherical uniformity
% 
%     F = gineajne(U)
%
%     A weighted sum of Gine's and Anje's statistics is consistent against
%     all alternative to uniformity on S^(p-1), the unit sphere in R^p.
%
%     INPUTS
%     U - [n x p] matrix, n samples with dimensionality p
%         the data should already be projected to the unit hypersphere
%
%     OUTPUTS
%     F - statistic
%
%     REFERENCE
%     Prentice, MJ (1978). On invariant tests of uniformity for directions
%       and orientations. Annals of Statistics 6: 169-176.
%
%     SEE ALSO
%     spheresign

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     The full license and most recent version of the code can be found on GitHub:
%     https://github.com/brian-lau/spheretest
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function F = gineajne(U)

[n,p] = size(U);

psi = sphere.psivec(U,n);
G = n/2 - (p-1)/(2*n) * (gamma((p-1)/2)/gamma(p/2))^2 * sum(sin(psi));
A = (n/4) - (1/(n*pi))*sum(psi);
F = G + A;
