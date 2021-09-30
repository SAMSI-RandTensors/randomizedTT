function X = TTrounding(Y, tol, max_rank)
% TTROUNDING  rounding procedure for TT-tensors with right to left orthogonalization followed by left to right compression steps.
%
%  X = TTROUNDING(Y, tol) rounds the TT-tensor Y with tolerance tol, yielding a new TT-tensor X with left-orthogonal cores X{1}, .., X{N-1} and smaller ranks.
%
%  X = TTROUNDING(Y, tol, max_rank) rounds the TT-tensor Y with tolerance tol and maximum bound on TT-ranks max_rank, yielding a new TT-tensor X with left-orthogonal cores X{1}, .., X{N-1} and smaller ranks.

if (nargin == 2)
    max_rank = 0;
end

X = TTorthogonalizeRL(Y);
[N,I,~] = TTsizes(X);

normX = norm(X{1}, 'fro');
tau = tol * normX / sqrt(N - 1);    % Truncation threshold

for n = 1:N-1
    [X{n},R] = qr(X{n}, 0);
    [U, S, V] = svd(R, 'econ');
    rank = trunc(diag(S), tau, max_rank);
    X{n} = X{n}*U(:, 1:rank);
    X{n+1} = S(1:rank, 1:rank) * V(:, 1:rank)' * v2h(X{n+1}, I(n+1));
    X{n+1} = h2v(X{n+1}, I(n+1));
end