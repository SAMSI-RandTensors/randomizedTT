function X = TTKronOp(A, Y)
% KRONOP  Apply an operator in Kronecker form to a TT-tensor.
%
% X = KronOp(A, Y)  returns the result of applying an operator in Kronecker form A{1} ⊗ ... ⊗ A{N} to a TT-tensor X. 

[N,I,ry] = TTsizes(Y);

X = cell(N, 1);
for n = 1:N
    assert(I(n) == size(A{n},2), "incompatible mode dimensions in TTKronOp")
    
    y = reshape(Y{n}, ry(n), I(n), ry(n+1));
    x = zeros(ry(n), size(A{n},1), ry(n+1));
    for r=1:ry(n)
        x(r,:,:) = A{n} * reshape(y(r,:,:), I(n), ry(n+1));
    end
    X{n} = reshape(x, [], ry(n+1));
end