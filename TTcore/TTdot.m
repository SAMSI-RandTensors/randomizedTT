function w = TTdot(X, Y)
% TTDOT  Compute the inner product between two TT-tensors.
% 
% w = TTdot(X,Y) returns the inner product w between two TT-tensors X and Y with compatible mode sizes using successive contractions in a recursive manner from left to right.


[N,I,~] = TTsizes(X);

% forming W_n = contractions between cores (1, ..., n) of X and Y 
w = X{1}' * Y{1};
for n = 2:N
    w = X{n}' * h2v(w * v2h(Y{n}, I(n)), I(n));
end
