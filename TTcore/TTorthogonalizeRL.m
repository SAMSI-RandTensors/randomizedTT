function X = TTorthogonalizeRL(Y)
% TTORTHOGONALIZERL  left to right orthogonalization of a TT-tensor.
%
%  X = TTORTHOGONALIZERL(Y) returns a TT-tensor X equivalent to Y with right-orthogonal cores X{2}, ..., X{N}.

X = Y;
[N,I,~] = TTsizes(X);
X = cell(N,1);

% Right to Left orthogonalization of Y
X{N} = Y{N};
for n = N:-1:2
    [Q, R] = qr(v2h(X{n}, I(n))', 0);
    X{n} = h2v(Q', I(n));
    X{n-1} = Y{n-1} * R';
end