function k = distkern(X,Y)

Yt = Y';
XX = sqrt(sum(X.*X,2));
YY = sqrt(sum(Yt.*Yt));
D = sqrt(utils.sqdist(X,Y));

k = 0.5 * (bsxfun(@plus,XX,YY) - D);
