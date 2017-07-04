function Y = fwht(X)

[n,m] = size(X);
n2 = nextpow2(n);

% Zero-pad to nextpow2
if n ~= 2^n2
   X = [X ; zeros(2^n2-n,m)];
end

try
   % Scaled to match Matlab fwht
   Y = utils.mexHadamard(X)/2^n2;
catch
   Y = fwht(X,2^n2,'hadamard');
end

