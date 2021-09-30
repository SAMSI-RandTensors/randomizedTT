function [N,I,R] = TTsizes(X)
% TTSIZES  compute dimension, mode and rank sizes of a TT-tensor.
%
% [N,I,R] = TTSIZES(X) return the dimension N, 
%                             the mode sizes I(1), ..., I(N)
%                             and ranks R(1), ..., R(N+1) of TT-tensor X.

N = size(X, 1);                % number of modes
R = zeros(N + 1, 1);           % TT ranks of X
R(1) = 1;                      % Left boundary condition
R(N + 1) = 1;                  % Right boundary condition

I = zeros(N, 1);               % modes sizes
I(1) = size(X{1}, 1);

for n = 2:N
    R(n) = size(X{n - 1}, 2);
    I(n) = size(X{n}, 1)/R(n);
end