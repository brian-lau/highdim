function stat = hsic_(K,L,m,biased)

% Song et al. eq 4, 5
if biased
   H = eye(m) - ones(m)/m;
   Kc = K*H;
   Lc = L*H;
   stat = trace(Kc*Lc)/(m-1)^2;
else
   K = utils.zerodiag(K);
   L = utils.zerodiag(L);
   l = ones(m,1);
   h = trace(K*L) + (l'*K*l*l'*L*l)/(m-1)/(m-2) - 2*(l'*K*L*l)/(m-2);
   stat = h/m/(m-3);
end