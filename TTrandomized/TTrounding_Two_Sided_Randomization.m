function X = TTrounding_Two_Sided_Randomization(Y, l)
% TTROUNDING_TWO_SIDED_RANDOMIZATION rounding procedure for TT-tensors with no orthogonalization and two-sided randomized compression steps.
%
% X = TTROUNDING_TWO_SIDED_RANDOMIZATION(Y, l) rounds the TT-tensor Y, yielding a new TT-tensor X with prescribed ranks l(1), ..., l(N+1).
%
% See also TTROUNDING_ORTHOGONALIZE_THEN_RANDOMIZE and TTROUNDING_RANDOMIZE_THEN_ORTHOGONALIZE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements Algorithm 3.3: TT-rounding, Two-Sided-Randomization. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,I,~] = TTsizes(Y);

% Heuristic choice of oversampling ranks
rho = [1; floor(2*l(2:N)); 1];

WL = TTpartialContractionsLR(TTrand(I, l), Y);
WR = TTpartialContractionsRL(Y, TTrand(I, rho));

X = cell(N,1);
[U,S,V] = svd(WL{1} * WR{1}, 'econ');
S = diag(diag(S).^(-0.5));
L = WR{1}*V*S;
R = S*U'*WL{1};

X{1} = Y{1} * L;

for n=2:N-1
   [U,S,V] = svd(WL{n} * WR{n}, 'econ');
   S = diag(diag(S).^(-0.5));
   L = WR{n}*V*S;
   X{n} = h2v(R * v2h(Y{n} * L, I(n)), I(n));
   
   R = S*U'*WL{n};
end

X{N} = h2v( R * v2h(Y{N}, I(N)), I(N));