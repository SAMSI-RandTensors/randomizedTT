function X = TTrand(I, R)
% TTRAND  compute a random TT-tensor with normalized independent Gaussian-distributed entries.
%
% X = TTrand(I,R)  returns a Random Gaussian TT-tensor with mode sizes I(1), ..., I(N) and TT ranks R(1), ..., R(N + 1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements Definition 3.1 of a Random Gaussian TT-Tensor. %
% Entries of X{i} are normally distributed with mean 0      %
% and variance 1/(I(i)*R(i)*R(i+1)).                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(I);
if( length(R) == 1)
    R = [1; R * ones(N - 1, 1); 1];
end

assert(length(R) == N+1, "Please provide conforming sizes of the modes sizes and the TT ranks");
assert(R(1) == 1 && R(N+1) == 1, "Boundary TT ranks must be 1");

X = cell(N, 1);
for n = 1 : N
    X{n} = randn(I(n) * R(n), R(n + 1));
    s = sqrt(I(n)*R(n)*R(n+1));             % variance normalization
    X{n} = X{n}/s;
end