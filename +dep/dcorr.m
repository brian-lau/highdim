function d = dcorr(x,y,modified)

if nargin < 3
   modified = false;
end

if modified
   [d,dvx,dvy] = dep.dcov(x,y,true);
   d = d/sqrt(dvx*dvy);
else
   [d,dvx,dvy] = dep.dcov(x,y,false);
   d = sqrt(d/sqrt(dvx*dvy));
end
