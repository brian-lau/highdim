% DEPTEST2                    Interface for two-sample (in)dependence tests
%
%     Given a sample X1,...,Xm from a p-dimensional multivariate distribution,
%     and a sample Y1,...,Xn from a q-dimensiona multivariate distribution,
%     test one of the hypotheses:
%
%     H0 : X and Y are drawn from the same distribution
%     
%     using the following tests,
%        'mmd' - Maximal Mean Discrepancy
%        'energy'
%        'ks'
%
%     H0 : X and Y are mutually independent
%
%     using the following tests,
%        'dcorr' - distance correlation (default)
%        'rv' - RV coefficient
%        'hsic' - Hilbert-Schmidt Independence Criterion
%
%     PROPERTIES
%     x       - [m x p] matrix, m samples with dimensionality p
%     x       - [n x q] matrix, n samples with dimensionality q
%     m       - # of x samples
%     p       - # of x dimensions
%     n       - # of y samples
%     q       - # of y dimensions
%     test    - string (see above, default = 'dcorr')
%     params  - 
%     alpha   - alpha level (default = 0.05)
%     stat    - corresponding statistic
%     pval    - p-value
%     h       - boolean, 1 indicates rejection of null at alpha
%     runtime - elapsed time for running test, in seconds
%
%     EXAMPLE
%
%     REFERENCE
%     Gretton et al (2008). A kernel statistical test of independence. NIPS
%     Szekely et al (2007). Measuring and testing independence by correlation 
%       of distances. Ann Statist 35: 2769-2794
%     Szekely & Rizzo (2013). The distance correlation t-test of independence 
%       in high dimension. J Multiv Analysis 117: 193-213
%
%     SEE ALSO
%     DepTest1, UniSphereTest

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     The full license and most recent version of the code can be found
%     at:help
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

classdef DepTest2 < hgsetget
   properties
      x
      y
   end
   properties (Dependent=true,SetAccess=private)
      m
      p
      n
      q
   end
   properties
      test
      params
      alpha = 0.05;
   end
   properties (SetAccess=private)
      stat
      pval
      h
      runtime
   end
   properties (Hidden=true,SetAccess=private)
      autoRun
      validTests = {'dcorr' 'rv' 'hsic'...
         'mmd'};
   end
   methods
      function self = DepTest2(varargin)
%          if (nargin == 1) || (rem(nargin,2) == 1)
%             varargin = {'x' varargin{:}};
%          end
         
         par = inputParser;
         par.KeepUnmatched = true;
         addParamValue(par,'x',[],@isnumeric);
         addParamValue(par,'y',[],@isnumeric);
         addParamValue(par,'autoRun',true,@islogical);
         addParamValue(par,'test','dcorr',@ischar);
         parse(par,varargin{:});

         self.autoRun = par.Results.autoRun;
         self.params = par.Unmatched;
         self.test = par.Results.test;
         self.replaceData(par.Results.x,par.Results.y);
      end
      
      function replaceData(self,x,y)
         old = self.autoRun;
         self.autoRun = false;
         self.x = x;
         self.y = y;
         self.autoRun = old;
         if ~isempty(self.x) && ~isempty(self.y) && self.autoRun
            self.run();
         end
      end
      
      function set.x(self,x)
         self.x = x;
         if ~isempty(self.x) && ~isempty(self.y) && self.autoRun
            self.run();
         end
      end
      
      function set.y(self,y)
         self.y = y;
         if ~isempty(self.x) && ~isempty(self.y) && self.autoRun
            self.run();
         end
      end
      
      function set.test(self,test)
         test = lower(test);
         if any(strcmp(test,self.validTests))
            self.test = test;
            if ~isempty(self.x) && ~isempty(self.y) && self.autoRun
               self.run();
            end
         else
            error('Invalid test');
         end
      end
      
      function set.params(self,params)
         self.params = params;
         if ~isempty(self.x) && self.autoRun
            self.run();
         end
      end
      
      function set.alpha(self,alpha)
         assert((alpha>0)&&(alpha<1),'0<alpha<1');
         self.alpha = alpha;
      end
            
      function m = get.m(self)
         m = size(self.x,1);
      end
      
      function n = get.n(self)
         n = size(self.y,1);
      end
      
      function p = get.p(self)
         p = size(self.x,2);
      end
      
      function q = get.q(self)
         q = size(self.y,2);
      end
      
      function h = get.h(self)
         h = self.pval<self.alpha;
      end
      
      function run(self)
         tic;
         switch self.test
            case {'dcorr'}
               [self.pval,self.stat] = ...
                  dep.dcorrtest(self.x,self.y,self.params);
            case {'hsic'}
               [self.pval,self.stat] = ...
                  dep.hsic(self.x,self.y,self.params);
            case {'rv'}
               [self.pval,self.stat] = ...
                  dep.rvtest(self.x,self.y);
            case {'mmd'}
               [self.pval,self.stat] = ...
                  diff.mmdtest(self.x,self.y,self.params);
            otherwise
               % Never
         end
         self.runtime = toc;
      end
   end
end