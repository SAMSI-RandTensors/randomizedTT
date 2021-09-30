function X = v2h(Y, I) 
% V2H  Reshape core from vertical to horizontal unfolding.
%
%  X = V2H(Y, I) returns the horizontal unfolding of core Y 
%                with mode size I.
%
% See also H2V.

X = reshape(Y, [size(Y, 1)/I, I * size(Y, 2)]);