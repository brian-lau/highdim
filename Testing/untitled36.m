clear all
s = 5;
p = 10;
n = 100;
d = [6 3 1 0.4 0.2];

same = randn(p,s);
V1 = [same , randn(p,p-s)];
V1 = utils.mgs(V1);
V2 = [same , randn(p,p-s)];
V2 = utils.mgs(V2);

D1 = diag([d,2*rand(1,p-s)]);
D2 = diag([d,2*rand(1,p-s)]);

S1 = V1*D1*D1*V1';
S2 = V2*D2*D2*V2';

x = randn(n,p)*V1*D1;
y = randn(n,p)*V2*D2;


for i = 1:100
   x = randn(n,p);
   y = randn(n,p);
   [aInd, pval, varToThisDimInTSpace, varToThisDimInSSpace, varToThisDimWorstCase, randCI, sigs] = ...
      varAlignment(x,y,5);
   sig(i,:) = sigs;
   pv(i) = pval;
end

