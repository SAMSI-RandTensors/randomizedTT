function X = TTscale(Y, alpha)
% TTSCALE  scale a tensor in TT-format.
%
%  Y = TTSCALE(X, alpha) returns a new TT-tensor Y representing the tensor X multiplied entry-wise by the scalar alpha.



X = Y;
X{1} = X{1} * alpha;