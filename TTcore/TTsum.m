function X = TTsum(S, a, tol)
% TTSUM  compute and rounds a linear combination of TT-tensors.
%
%  X = TTSUM(S, a) rounds the linear combination a(1) S{1} + ... + a(m) S{m} yielding a new TT-tensor X with right-orthogonal cores X{2}, ..., X{N}.
%
%  X = TTSUM(a, S, tol) rounds the linear combination a(1) S(1) + ... + a(M) S(M) with tolerance tol in Frobenius norm, yielding a new TT-tensor X with right-orthogonal cores X{2}, ..., X{N}.

if (nargin == 2)
    tol = 1e-10;
end

X = TTscale(S{1}, a(1));
for j=2:length(S)
   X = TTaxby(1.,X,a(j),S{j}); 
end
X = TTrounding(X, tol);
