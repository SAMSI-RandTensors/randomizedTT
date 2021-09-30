function X = TTrounding_Orthogonalize_then_Randomize(Y, l)
% TTROUNDING_ORTHOGONALIZE_THEN_RANDOMIZE rounding procedure for TT-tensors with right to left orthonormnalization and left to right randomized compression steps.
%
% X = TTROUNDING_ORTHOGONALIZE_THEN_RANDOMIZE(Y, l) rounds the TT-tensor Y, yielding a new TT-tensor X with prescribed ranks l(1), ..., l(N+1).
%
% See also TTROUNDING_RANDOMIZE_THEN_ORTHOGONALIZE and TTROUNDING_TWO_SIDED_RANDOMIZATION.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements Algorithm 3.1: TT-rounding, Orthogonalize-then-Randomize. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = TTorthogonalizeRL(Y);
[N,I,rx] = TTsizes(X);

% Rounding the tensor cores form left to right
for n = 1 : N-1
    Z = X{n};
    [X{n},~] = qr(Z * randn(rx(n+1), l(n+1)), 0);
    M = X{n}'*Z;
    X{n+1} = h2v(M * v2h(X{n+1}, I(n+1)), I(n+1));
end