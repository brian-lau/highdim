% RVTEST                      Test RV coefficient of dependence
% 
%     [pval,rv,nrv] = rvtest(x,y)
%
%     INPUTS
%     x - [n x p] n samples of dimensionality p
%     y - [n x q] n samples of dimensionality q
%
%     OUTPUTS
%     pval - p-value from Pearson type III approximation
%     rv   - RV coefficient
%     nrv  - normalized RV coefficient
%
%     REFERENCE
%     Josse et al (2008). Testing the significance of the RV coefficient.
%       Computational Statistics and Data Analysis 53: 82-91
%
%     SEE ALSO
%     rv, dcorr, dcorrtest

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     The full license and most recent version of the code can be found at:
%     https://github.com/brian-lau/highdim
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

function [pval,rv,nrv] = rvtest(x,y)

[n,~] = size(x);
if n ~= size(y,1)
   error('RVTEST requires x and y to have the same # of samples');
end

[rv,xx,yy] = dep.rv(x,y);

% mean
bx = trace(xx)^2/trace(xx^2);
by = trace(yy)^2/trace(yy^2);
mu_rv = sqrt(bx*by)/(n-1);
% variance
tx = (n-1)/((n-3)*(n-1-bx)) * ...
     (n*(n+1)*(sum(diag(xx).^2)/trace(xx^2)) - (n-1)*(bx+2));
ty = (n-1)/((n-3)*(n-1-by)) * ...
     (n*(n+1)*(sum(diag(yy).^2)/trace(yy^2)) - (n-1)*(by+2));
var_rv = (2*(n-1-bx)*(n-1-by))/((n+1)*(n-1)^2*(n-2)) *...
     (1 + ((n-3)/(2*n*(n-1)))*tx*ty);

nrv = (rv - mu_rv)/sqrt(var_rv);

a = asym;
if a >= 0
   pval = 1 - gamcdf(nrv-(-2/a),4/a^2,a/2);
else
   pval = gamcdf(a/abs(a)*nrv+2/abs(a),4/a^2,abs(a)/2);
end

% Skewness estimate for Pearson III approximation
% Translated from coeffRV.R in FactoMineR 1.2.7
function a = asym
   S3 = sum(diag(xx).^3);
   S3etoile = sum(diag(yy).^3);  
   U = sum(sum(xx.^3));
   Uetoile = sum(sum(yy.^3));  
   B = diag(xx)'*xx*diag(xx);
   Betoile = diag(yy)'*yy*diag(yy);  
   R = diag(xx)'*diag(xx*xx);
   Retoile = diag(yy)'*diag(yy*yy); 
   TT = sum(diag(xx));
   Tetoile = sum(diag(yy)); 
   S2 = sum(diag(xx).^2);
   S2etoile = sum(diag(yy).^2); 
   T2 = sum(diag(xx*xx));
   T2etoile = sum(diag(yy*yy));
   T3 = sum(diag(xx*xx*xx));
   T3etoile = sum(diag(yy*yy*yy));
   
   total = n^2*(n+1)*(n^2+15*n-4)*S3*S3etoile+...
      4*(n^4-8*n^3+19*n^2-4*n-16)*U *Uetoile+...
      24*(n^2-n-4)*(U*Betoile + B*Uetoile)+...
      6*(n^4-8*n^3+21*n^2-6*n-24)* B*Betoile+...
      12*(n^4-n^3-8*n^2+36*n-48)* R*Retoile+...
      12*(n^3-2*n^2+9*n-12)  *(TT*S2*Retoile+R*Tetoile*S2etoile) +...
      3*(n^4-4*n^3-2*n^2+9*n-12)*(TT*Tetoile*S2*S2etoile )  +...
      24*( (n^3-3*n^2-2*n+8)*(R*Uetoile+U*Retoile)+(n^3-2*n^2-3*n+12)*...
      (R*Betoile+B*Retoile))+...
      12*(n^2-n+4)*(TT*S2*Uetoile+U*Tetoile*S2etoile)+...
      6*(2*n^3-7*n^2-3*n+12)*(TT*S2*Betoile+B*Tetoile*S2etoile) -... 
      2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3etoile+(2*Uetoile+3*Betoile)*S3)-...
      3*n*((n-1)^2) * (n+4)*((TT*S2+4*R)*S3etoile+(Tetoile*S2etoile+4*Retoile)*S3)+...
      2*n*(n-1)*(n-2)*( (TT^3+6*TT*T2+8*T3)*S3etoile +...
      (Tetoile^3+6*Tetoile*T2etoile+8*T3etoile)*S3)  +...
      TT^3*((n^3-9*n^2+23*n-14)*Tetoile^3+ 6*(n-4)*Tetoile*T2etoile+8*T3etoile)+...
      6*TT*T2*((n-4)*Tetoile^3+(n^3-9*n^2+24*n-14)*Tetoile*T2etoile+4*(n-3)*T3etoile)+...
      8*T3*(Tetoile^3+3*(n-3)*Tetoile*T2etoile+(n^3-9*n^2+26*n-22)*T3etoile)  -...
      16*(TT^3*Uetoile+U*Tetoile^3)-6*(TT*T2*Uetoile+U*Tetoile*T2etoile)*(2*n^2-10*n+16)-...
      8*(T3*Uetoile+U*T3etoile)*(3*n^2-15*n+16)-(TT^3*Betoile+B*Tetoile^3)*(6*n^2-30*n+24)-...
      6*(TT*T2*Betoile+B*Tetoile*T2etoile)*(4*n^2-20*n+24)-...
      8*(T3*Betoile+B* T3etoile)*(3*n^2-15*n+24)   -...
      (n-2)*(24*(TT^3*Retoile+R*Tetoile^3)+6*(TT*T2*Retoile+R*Tetoile*T2etoile)*(2*n^2-10*n+24)+...
      8*(T3*Retoile+R*T3etoile)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(TT^3*Tetoile*S2etoile+TT*S2*Tetoile^3)+...
      6*(TT*T2*Tetoile*S2etoile+TT*S2*Tetoile*T2etoile)*(n^2-5*n+6)+...
      48*(T3*Tetoile*S2etoile+TT*S2*T3etoile));

   esperancet3 = total/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5));
   esperance = TT*Tetoile/(n-1);
   variance = (2*((n-1)*T2-TT^2)*((n-1)*T2etoile-Tetoile^2)/...
      (((n-1)^2)*(n+1)*(n-2))) + (n*(n+1)*S2-(n-1)*(TT^2+2*T2))*...
      (n*(n+1)*S2etoile-(n-1)*(Tetoile^2+2*T2etoile))/((n+1)*n*(n-1)*(n-2)*(n-3));
   cumulant3 = esperancet3-3*esperance *variance-esperance^3;
   a = cumulant3/(variance^(3/2));
end

end