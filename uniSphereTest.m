% rayleigh - not consistent against alternatives with zero resultant length
%            Mardia & Jupp, pg 209
% most powerful invariant test against von mises alternative
% bingham  - antipodially symmetric
%            not consistent alterriatives with E[xx?] = (1/p)*Ip
function [pval,stat] = uniSphereTest(x,varargin)

if nargin == 2
   [pval,stat] = uniSphereTest(x,'test',varargin{1});
   return
elseif rem(nargin,2) == 0
   [pval,stat] = uniSphereTest(x,'test',varargin{1},varargin{2:end});
   return
end

par = inputParser;
addRequired(par,'x',@isnumeric);
addParamValue(par,'test','rayleigh',@ischar);
addParamValue(par,'nboot',1000,@isnumeric);
addParamValue(par,'k',25,@isnumeric);
addParamValue(par,'dist','boot',@ischar);
parse(par,x,varargin{:});

[n,p] = size(x);
U = spheresign(x);

switch lower(par.Results.test)
   case {'rayleigh','r'}
      stat = rayleigh(U,n,p);
      pval = 1 - chi2cdf(stat,p);
   case {'gine','g'}
      stat = gine(U,n,p);
      [pval,boot] = bootstrap('gine',stat,par.Results.nboot,n,p);
   case {'gine3','g3'}
      if p ~= 3
         error('Only valid for p = 3');
      end
      stat = gine3(U,n,p);
      pval = 1 - sumchi2cdf(stat,3);
   case {'ajne','a'}
      stat = ajne(U,n);
      [pval,boot] = bootstrap('ajne',stat,par.Results.nboot,n,p);
   case {'gine-ajne','ga','ag'}
      stat = gineajne(U,n,p);
      [pval,boot] = bootstrap('gineajne',stat,par.Results.nboot,n,p);
   case {'bingham','b'}
      stat = bingham(U,n,p);
      pval = 1 - chi2cdf(stat,((p-1)*(p+2))/2);
   case {'rp'}
      pval = rp(U,n,p,par.Results.k,...
         par.Results.dist,par.Results.nboot);
      [~,~,adj_p] = fdr_bh(pval,.05,'pdep');
      stat = [];%pval;
      pval = min(adj_p);
   otherwise
      error('Unknown test.');
end

function [pval,boot] = bootstrap(f,stat,nboot,n,p)
boot = zeros(nboot,1);
for i = 1:nboot
   Umc = spheresign(randn(n,p));
   boot(i) = feval(f,Umc,n,p);
end
pval = sum(boot>=stat)/nboot;

function Rnp = rayleigh(U,n,p)
if 0
   Rnp = (p/n)*sum(sum(U*U'));
else
   % Modified Rayleigh test statistic (Mardia & Jupp, 10.4)
   Ubar = mean(U);
   T = n*p*sum(Ubar.^2);
   Rnp = (1-1/(2*n))*T + (1/(2*n*(p+2)))*T^2;
end

function psi = psivec(U,n)
xx = triu(U*U',1);
ind = triu(ones(n,n),1);
psi = acos(xx(ind==1));

function G = gine(U,n,p)
psi = psivec(U,n);
G = n/2 - (p-1)/(2*n) * (gamma((p-1)/2)/gamma(p/2))^2 *...
      sum(sin(psi));

function Fn = gine3(U,n,p)
psi = psivec(U,n);
Fn = (3*n)/2 - (4/(n*pi)) * sum(psi + sin(psi));

function A = ajne(U,n,p)
psi = psivec(U,n);
A = (n/4) - (1/(n*pi))*sum(psi);

function F = gineajne(U,n,p)
psi = psivec(U,n);
G = n/2 - (p-1)/(2*n) * (gamma((p-1)/2)/gamma(p/2))^2 *...
      sum(sin(psi));
A = (n/4) - (1/(n*pi))*sum(psi);
F = G + A;

function B = bingham(U,n,p)
if 1
   T = (1/n)*U'*U;
   B = ((n*p*(p+2))/2)*(trace(T^2) - 1/p);
else
   % Modified Bingham test statistic (Mardia & Jupp, 10.7)
   % seems to blow up for certain data?
   T = (1/n)*U'*U;
   B = ((n*p*(p+2))/2)*(trace(T^2) - 1/p);
   B0 = (2*p^2+3*p+4)/(6*(p+4));
   B1 = -(4*p^2+3*p-4)/(3*(p+4)*(p^2+p+2));
   B2 = 4*(p^2-4)/(3*(p+4)*(p^2+p+2)*(p^2+p+6));
   B = B*(1 - (1/n)*(B0 + B1*B + B2*B^2));
end

function pval = rp(U,n,p,k,dist,nboot)
u0 = spheresign(randn(k,p));
Y = zeros(n,k);
pval = zeros(k,1);
switch dist
   case 'asymp'
      for i = 1:k
         Y(:,i) = acos(U*u0(i,:)');
         test_cdf = [ Y(:,i) , rpcdf(Y(:,i),p)];
         [~,pval(i)] = kstest(Y(:,i),'CDF',test_cdf);
      end
   otherwise
      for i = 1:k
         Y(:,i) = U*u0(i,:)';
      end
      Umc = spheresign(randn(nboot,p));
      Ymc = Umc*u0(1,:)';
      for i = 1:k
         [~,pval(i)] = kstest2(Y(:,i),Ymc);
      end
end


