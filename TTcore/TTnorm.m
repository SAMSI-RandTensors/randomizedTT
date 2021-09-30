function normx = TTnorm(X, method)
% TTNORM  compute the Frobenius norm of a TT-tensor.
%
%  normx = TTNORM(X) computes the Frobenius norm of TT-tensor X using the inner product strategy (Gram matrices)
%
%  normx = TTNORM(X, "OLR") computes the Frobenius norm of TT-tensor X using the left to right orthogonalization strategy
%
%  normx = TTNORM(X, "ORL") computes the Frobenius norm of TT-tensor X using the right to left orthogonalization strategy
%
%  normx = TTNORM(X, "G") is the same as normx = TTNORM(X, "G")

if(nargin > 1)
    method = method{1};
    if(~strcmp(method, "OLR") && ~strcmp(method, "ORL") && ~strcmp(method, "G"))
        fprintf("unknown norm computation strategy\n");
    end
else
    method = 'G';
end

[N,I,~] = TTsizes(X);

if(strcmp(method, "OLR"))
    % Left to right orthogonalization of Y
    for n = 1 : N - 1
        [X{n}, R] = qr(X{n}, 0);
        X{n + 1} = h2v(R * v2h(X{n + 1}, I(n + 1)), I(n + 1));
    end
    
    normx = norm(X{N}, 'fro');
elseif(strcmp(method, "ORL"))
    % Right to left orthogonalization of Y
    for n = N : -1 : 2
        [X{n}, R] = qr(v2h(X{n}, I(n))', 0);
        X{n} = h2v(X{n}', I(n));
        X{n - 1} = X{n - 1} * R';
    end
    
    normx = norm(X{1}, 'fro');
else
    G = X{1}' * X{1};
    for n = 2 : N
        G = X{n}' * h2v(G * v2h(X{n}, I(n)), I(n));
        G = (G + G')/2;
    end
    normx = sqrt(abs(G));
end