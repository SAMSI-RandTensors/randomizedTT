function X = h2v(Y, I) 
% H2V  Reshape core from vertical to horizontal unfolding.
%
%  X = H2V(Y, I) returns the horizontal unfolding of core Y with mode size I.
%
% See also V2H.

X = reshape(Y, [size(Y, 1) * I, size(Y, 2)/I]);