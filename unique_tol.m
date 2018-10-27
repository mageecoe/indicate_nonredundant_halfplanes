function [Amin, ind_nred] = unique_tol(A, tol)
%UNIQUE_TOL Like uniquetol, with well-defined rows to keep
%
% [Amin, IND_NRED] = UNIQUE_TOL(A, TOL) removes duplicate rows of A,
% keeping the last occuring row of the duplicates. TOL is a numerical
% tolerance used for determining rows as equal. IND_NRED is a logical
% array whose true indices represent the set of unique rows of A.
%
% % Example
% % -------
%
% A = [1 0 0 0;...
%      0 2 0 0;...
%      0 0 3 0;...
%      1 0 0 0;...
%      0 0 0 4;...
%      0 2 0 0];
% [Amin, ind_nred] = unique_tol(A);
%
% assert(all(all(Amin == A([3 4 5 6],:))))

if nargin < 2
    tol = 1e-5;
end

if not(isscalar(tol))
    error('The second argument, tol, must be a scalar');
end

% Initialize
[m,n] = size(A);
ind_nred = true(m,1);
rows_to_check = true(m,1);

% Identify duplicate rows
max_it = m;
for it=1:max_it


    % Find row
    if any(rows_to_check)
        j = find(rows_to_check,1);
        rows_to_check(j) = false;
    else
        break
    end

    % Indices of possible duplicate rows
    ind_duplicate_candidates = rows_to_check;

    % Loop over columns
    for i=1:n

        % Check stopping criterion
        if not(any(ind_duplicate_candidates))
            break
        end

        A_col = A(:,i) - A(j,i);

        %n_too_large = find(abs(A_col) > tol);
        %ind_duplicate_candidates(n_too_large) = false;
        
        % The following rows give the same result as above, but is
        % significantly faster
        ind_too_large = abs(A_col) > tol;
        ind_duplicate_candidates = and(ind_duplicate_candidates, ~ind_too_large);
        
    end

    if any(ind_duplicate_candidates)
        ind_nred(j) = false;
    end

end

Amin = A(ind_nred,:);

% end unique_tol
