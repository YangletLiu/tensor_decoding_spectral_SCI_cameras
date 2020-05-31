function X = TRNNM(Y,tau,R)

[U,S,V] = svd(Y,'econ');
s = diag(S);
W = 1./(abs(s) + eps);
W(1:R) = 0;
s = diag(max(s - tau.*W,zeros(size(s))));
X = U*s*V';
