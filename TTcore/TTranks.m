function R = TTranks(X)
% TTRANKS  compute the ranks of a TT-tensor.
%
% R = TTRANKS(X) returns the TT-ranks R(1), ..., R(N+1) of TT-tensor X. 

N = length(X);
R = zeros(N + 1, 1);
R(1) = 1;
for n = 1:N
    R(n+1) = size(X{n}, 2);
end