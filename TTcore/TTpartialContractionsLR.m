function W = TTpartialContractionsLR(X, Y)
% TTPARTIALCONTRACTIONSLR  Compute the sequence of nested contractions between two TT-tensors from left to right.
%
%  W = TTPARTIALCONTRACTIONSLR(X, Y) returns matrices W{1}, ..., W{N-1} where W{i} is the contraction between cores (1, ..., i) of X and Y.
%
% See also TTPARTIALCONTRACTIONSRL.

[N,I,~] = TTsizes(X);
W = cell(N-1, 1);

W{1} = X{1}' * Y{1};
for n = 2:N-1
    W{n} = X{n}' * h2v(W{n-1} * v2h(Y{n}, I(n)), I(n));
end

