function [pval,stat] = cpctest(x,s)

nS = numel(x);
for i = 1:nS
   n(i) = size(x{i},1);
   S(:,:,i) = cov(x{i});
end

      nboot = 500;
method = 'fuj'
switch method
   case {'boot'}
      [Q,D] = dim.cpca(S,n,'tol',1e-10);
      stat = lrt(D,S,n);
      
      z = cat(1,x{:});
      cumn = cumsum(n);
      for i = 1:nboot
         ind = randperm(sum(n));
         for j = 1:nS
            if j == 1
               start = 0;
            else
               start = cumn(j-1);
            end
            ind2 = (start+1):cumn(j);
            xb{j} = z(ind(ind2),:);
            Sb(:,:,j) = cov(xb{j});
         end
         [~,Db] = dim.cpca(Sb,n);
         t(i) = lrt(Db,Sb,n);
      end
      pval = sum(t>=stat)/nboot;
   case 'schott'
      s = 2;
      stat = schottPartial(S,n,s);
      
      z = cat(1,x{:});
      cumn = cumsum(n);
      for i = 1:nboot
         ind = randperm(sum(n));
         for j = 1:nS
            if j == 1
               start = 0;
            else
               start = cumn(j-1);
            end
            ind2 = (start+1):cumn(j);
            xb{j} = z(ind(ind2),:);
            Sb(:,:,j) = cov(xb{j});
         end
         ts(i) = schottPartial(Sb,n,s);
      end
      pval = sum(ts>=stat)/nboot;
   case {'fuj'}
      %s = 2;
      %[~,stat] = dim.kryzsim(x{1},x{2},s);
      [~,stat] = dim.krzsim(x{1},x{2},s);
      %stat = dim.kryzsim(x{1},x{2},s);
      
%       z = cat(1,x{:});
%       cumn = cumsum(n);
%       for i = 1:nboot
%          ind = randperm(sum(n));
%          for j = 1:nS
%             if j == 1
%                start = 0;
%             else
%                start = cumn(j-1);
%             end
%             ind2 = (start+1):cumn(j);
%             xb{j} = z(ind(ind2),:);
%          end
%          %[~,Tm(i)] = dim.kryzsim(xb{1},xb{2},s);
%          [~,Tm(:,i)] = dim.krzsim(xb{1},xb{2},s);
%       end
      %z = cat(2,x{:});
      p = size(x{1},2);
      for i = 1:nboot
         [~,Tm(:,i)] = dim.krzsim(randn(n(1),p),randn(n(2),p),s);
%          ind1 = reshape(randperm(n(1)*p),n(1),p);
%          ind2 = reshape(randperm(n(2)*p),n(2),p);
%          [~,Tm(:,i)] = dim.krzsim(x{1}(ind1),x{2}(ind2),s);
      end
      %pval = sum(Tm>=stat)/nboot
      for i = 1:numel(s)
         pval(i) = sum(Tm(i,:)<=stat(i))/nboot;
      end
   otherwise
end

function ts = schottPartial(S,n,s)
nS = size(S,3);
p = size(S,1);
P = zeros(p);
for i = 1:nS
   [Q,D] = eig(S(:,:,i));
   P = P + Q(:,1:s)*Q(:,1:s)';
end
[Q,D] = eig(P);
rho = diag(D);
ts = sum(n)*sum(rho(s+1:p));

function t = lrt(Dr,S,n)
k = size(S,3);
t = 0;
for i = 1:k
   D = diag(Dr(:,i));
   t = t + n(i)*log(det(D)/det(S(:,:,i)));
end