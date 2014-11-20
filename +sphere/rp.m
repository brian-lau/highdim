% RP                          Random projection statistic for spherical uniformity 
% 
%     stat = rp(U,k)
%
%     INPUTS
%     U - [n x p] matrix, n samples with dimensionality p
%         the data should already be projected to the unit hypersphere
%     k - number of random vectors to project onto
%
%     OUTPUTS
%     stat - [n x k] vector of of angles between data and k random vectors
%
%     REFERENCE
%     Cuesta-Albertos, JA et al (2009). On projection-based tests for directional 
%       and compositional data. Stat Comput 19: 367-380
%     Cuesta-Albertos, JA et al (2007). A sharp form of the Cramer?Wold theorem. 
%       J Theor Probab 20: 201?209
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

function stat = rp(U,k)

[n,p] = size(U);

% Uniform random directions
u0 = spheresign(randn(k,p));
stat = zeros(n,k);

for i = 1:k
   stat(:,i) = acos(U*u0(i,:)');
end
