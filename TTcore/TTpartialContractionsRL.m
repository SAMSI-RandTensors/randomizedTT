function W = TTpartialContractionsRL(X, Y)
% TTPARTIALCONTRACTIONSRL  Compute the sequence of nested contractions between two TT-tensors from right to left.
%
%  W = TTPARTIALCONTRACTIONSRL(X, Y) returns matrices W{1}, ..., W{N-1} where W{i} is the contraction between cores (i+1, ..., N) of X and Y.
%
% See also TTPARTIALCONTRACTIONSLR.

[N,I,~] = TTsizes(X);
W = cell(N-1, 1);

W{N-1} = v2h(X{N}, I(N)) * v2h(Y{N}, I(N))';
for n = N-1:-1:2
    W{n-1} = v2h(X{n} * W{n}, I(n)) * v2h(Y{n}, I(n))';
end

