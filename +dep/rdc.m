% Based on R code:
% https://github.com/lopezpaz/randomized_dependence_coefficient/blob/master/code/algorithms.r
% rdc <- function(x,y,k=20,s=1/6,f=sin) {
%   x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)),1)
%   y <- cbind(apply(as.matrix(y),2,function(u)rank(u)/length(u)),1)
%   x <- s/ncol(x)*x%*%matrix(rnorm(ncol(x)*k),ncol(x))
%   y <- s/ncol(y)*y%*%matrix(rnorm(ncol(y)*k),ncol(y))
%   cancor(cbind(f(x),1),cbind(f(y),1))$cor[1]
% }

function r = rdc(x,y,varargin)

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
addParamValue(par,'k',1.5,@isscalar);
addParamValue(par,'s',[],@isscalar);
addParamValue(par,'f',@sin,@(x) isa(x,'function_handle'));
addParamValue(par,'demean',false,@islogical);
parse(par,x,y,varargin{:});

n = size(x,1);
if par.Results.demean
   x = bsxfun(@minus,x,mean(x));
   y = bsxfun(@minus,y,mean(y));
end

x = [tiedrank(x)/n ones(n,1)];
y = [tiedrank(y)/n ones(n,1)];

x = par.Results.f(s/size(x,2)*x*randn(size(x,2),k));
y = par.Results.f(s/size(y,2)*y*randn(size(y,2),k));

warning off;
[~,~,r] = canoncorr([x ones(n,1)],[y ones(n,1)]);
warning on;
r = r(1);

