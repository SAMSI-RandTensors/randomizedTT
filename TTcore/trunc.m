function r = trunc(svals, tol, max_rank)
% TRUNC  compute the smallest rank given a list of singular values.
%
% r = TRUNC(svals, tol)  returns the smallest rank r achieving tolerance tol, measured in the Frobenius norm, given the ordered list of singular values svals.
%
% r = TRUNC(svals, tol, max_rank)  returns the smallest rank r achieving tolerance tol, measured in the Frobenius norm, and/or smaller than max_rank, given the ordered list of singular values svals.

if nargin == 2 || max_rank == 0
    max_rank = length(svals);
else
    max_rank = min(length(svals), max_rank);
end


ss = sqrt(cumsum(svals.^2, 'reverse'));

for j = 2:max_rank
    if(ss(j) < tol)
        r = j - 1;
        return
    end
end

r = max_rank;