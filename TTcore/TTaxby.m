function X = TTaxby(alpha, Y, beta, Z)
% TTAXBY  compute the linear combination of two TT-tensors.
%
%  X = TTAXBY(alpha, Y, beta, Z) returns the formal TT-representation X of the linear combination alpha Y + beta Z with compatible TT-tensors Y, Z and coefficients alpha, beta.

[N,I,ry] = TTsizes(Y);
[~,~,rz] = TTsizes(Z);

rx = ry + rz;
rx(1) = 1;
rx(N+1) = 1;

X = cell(N, 1);         % initiate the result
X{1} = [alpha * Y{1}, beta * Z{1}];

for n = 2 : N-1
    X{n} = zeros(rx(n), I(n), rx(n+1));
    
    Yn = reshape(Y{n}, ry(n), I(n), ry(n+1));
    Zn = reshape(Z{n}, rz(n), I(n), rz(n+1));
    
    X{n}(1:ry(n), :, 1:ry(n+1)) = Yn;
    X{n}(ry(n)+1:ry(n)+rz(n), :, ry(n+1)+1:ry(n+1)+rz(n+1)) = Zn; 
    
    X{n} = reshape(X{n}, rx(n)*I(n), rx(n+1));
end

YN = reshape(Y{N}, ry(N), I(N));
ZN = reshape(Z{N}, rz(N), I(N));
X{N} = reshape([YN; ZN], rx(N)*I(N), rx(N+1));
