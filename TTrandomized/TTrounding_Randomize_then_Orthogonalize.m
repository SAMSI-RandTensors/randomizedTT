function X = TTrounding_Randomize_then_Orthogonalize(Y, l)
% TTROUNDING_RANDOMIZE_THEN_ORTHOGONALIZE rounding procedure for TT-tensors with no orthogonalization and left to right randomized compression steps.
%
% X = TTROUNDING_RANDOMIZE_THEN_ORTHOGONALIZE(Y, l) rounds the TT-tensor Y, yielding a new TT-tensor X with prescribed ranks l(1), ..., l(N+1) and left-orthogonal cores X{1}, ..., X{N-1}.
%
% See also TTROUNDING_ORTHOGONALIZE_THEN_RANDOMIZE and TTROUNDING_TWO_SIDED_RANDOMIZATION.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements Algorithm 3.2: TT-rounding, Randomize-then-Orthogonalize. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,I,~] = TTsizes(Y);
W = TTpartialContractionsRL(Y, TTrand(I, l));

X = cell(N,1);
X{1} = Y{1};
for n = 1 : N - 1
    Z = X{n};
    [X{n}, ~] = qr(Z * W{n}, 0);
    M = X{n}' * Z;
    X{n + 1} = h2v(M * v2h(Y{n + 1}, I(n+1)), I(n+1));
end
