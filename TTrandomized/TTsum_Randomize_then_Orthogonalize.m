function X = TTsum_Randomize_then_Orthogonalize(S, a, tol)
% TTSUM_RANDOMIZE_THEN_ORTHOGONALIZE  randomized summation and rounding procedure for a linear combination of TT-tensors.
%
% X = TTSUM_RANDOMIZE_THEN_ORTHOGONALIZE(S, a) rounds the linear combination a(1) S{1} + ... + a(m) S{m} yielding a new TT-tensor X with right-orthogonal cores X{2}, ..., X{N}.
%
% X = TTSUM_RANDOMIZE_THEN_ORTHOGONALIZE(S, a, tol) rounds the linear combination a(1) S{1} + ... + a(m) S{m} with tolerance tol in Frobenius norm, yielding a new TT-tensor X with right-orthogonal cores X{2}, ..., X{N}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements Algorithm 3.4: TT-rounding of a Sum, Randomize-then-Orthogonalize. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 2)
    tol = 1e-10;
end

m = length(S);
assert(length(a)==m, "Incompatible size of coefficients and summands in TTrounding_Sum");

[N,I,r1] = TTsizes(S{1});
rs = zeros(N+1,m);
rs(:,1) = r1;
for j=2:m
   rs(:,j) = TTranks(S{j}); 
end

% heuristic for the target ranks for the randomization procedure
ro = max(rs, [], 2) * 2 + m;
ro(1) = 1;
ro(N + 1) = 1;

O = TTrand(I, ro);

% Compute partial random contractions of summands from right to left
W = cell(N-1,m);
for j = 1:m
    S{j} = TTscale(S{j}, a(j));
    W(:,j) = TTpartialContractionsRL(S{j}, O);
end

% Form the complete random projections
for n = 1:N-1
    W{n,1} = vertcat(W{n,:});
    for j=2:m
        W{n,j} = [];
    end
end
W = W(:,1);

% Initialize with the first core of the formal TT-sum
X = cell(N, 1);
lr = [0, cumsum(rs(2, :))];
X{1} = zeros(I(1), lr(m+1));
for j = 1:m
    X{1}(:, lr(j)+1 : lr(j+1)) = S{j}{1};
end

% Left-to-right randomization and orthogonalization
for n = 1:N-1
    Zn = X{n};
    Yn = X{n} * W{n};
    [X{n}, ~] = qr(Yn, 0);
    Mn = X{n}' * Zn;
    lr = [0, cumsum(rs(n + 1, :))];
    rr = [0, cumsum(rs(n + 2, :))];
    if(n < N - 1)
        X{n+1} = zeros(size(Mn,1) * I(n + 1), rr(m + 1));
        for j = 1:m
            x =  Mn(:, lr(j)+1:lr(j+1)) * v2h(S{j}{n+1}, I(n+1));
            X{n+1}(:, rr(j)+1:rr(j+1)) = h2v(x, I(n+1));
        end
    else
        X{N} = Mn(:, 1:lr(2)) * v2h(S{1}{N}, I(N));
        for j = 2:m
            X{N} = X{N} + Mn(:, lr(j)+1:lr(j+1)) * v2h(S{j}{N}, I(N));
        end
        X{N} = h2v(X{N}, I(N));
    end
end

% TT-rounding of Y, which is already orthogonalized.

normX = norm(X{N}, 'fro');
tau = tol * normX / sqrt(N - 1);                    % Truncation threshold
for n = N:-1:2
    [Q,R] = qr(v2h(X{n}, I(n))', 0);
    [U, S, V] = svd(R, 'econ');
    rank = trunc(diag(S), tau);
    X{n} = h2v( (Q*U(:, 1:rank))', I(n));
    X{n-1} = X{n-1} * (V(:, 1:rank) * S(1:rank, 1:rank));
end