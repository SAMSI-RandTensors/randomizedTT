function X = TTsummandsKronOp(A, Y)
% TTSUMMANDSKRONOP   form the summands of the application to a TT-tensor of an operator in Kronecker decomposition form.
%
% X = TTSUMMANDSKRONOP(A, Y)  computes an array of TT-tensors X{j} with cores
%                 X{j}{n} = A{n,j} x₂ Y{j}{n}
%          such that the application of the operator A is given by X = Σ X{j}.


rA = size(A,2);
X = cell(rA, 1);
for j = 1 : rA
    X{j} = TTKronOp(A(:,j), Y);
end