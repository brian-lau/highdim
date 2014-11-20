% RPPDF                       Distribution of angles on a uniform hypersphere
% 
%     h = rppdf(theta,p)
%
%     The distribution of pairwise angles between vectors X1, ? ? ? ,Xn 
%     that are random points independently chosen with the uniform distribution 
%     on S^(p?1), the unit sphere in R^p.
%
%     INPUTS
%     theta - angles (radians) to evaluate pdf
%     p     - dimensionality (R^p)
%
%     OUTPUTS
%     h     - pdf
%
%     REFERENCE
%     Cai, T et al (2013). Distribution of angles in random packing on
%     spheres. J of Machine Learning Research 14: 1837-1864.
%
%     SEE ALSO
%     rp

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

function h = rppdf(theta,p)

h = (1/sqrt(pi)) * (gamma(p/2)/(gamma((p-1)/2)))*...
   (sin(theta).^(p-2));
