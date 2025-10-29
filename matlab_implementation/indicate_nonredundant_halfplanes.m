function [Amin, bmin, ind_minrep, ind_not_verified] = indicate_nonredundant_halfplanes(A, b, ind_to_check, z)
%INDICATE_NONREDUNDANT_HALFPLANES  Minimal H-representation algorithm
%
% [Amin, bmin] = indicate_nonredundant_halfplanes(A, b) returns a minimal
% representation of the polytope {x | A x <= b}.
%
% For details about the algorithm, see the paper
%
% Emil Klintberg, Magnus Nilsson, Lars Johannesson Mårdh. "A primal
% active-set minimal representation algorithm for polytopes with
% application to invariant-set calculations".
%
% The general syntax of the function is
%
% [Amin, bmin, ...
%  IND_MINREP, IND_NOT_VERIFIED] = indicate_nonredundant_halfplanes(A, b, IND_TO_CHECK, z)
%
% The input z is an interior point of the polytope, i.e. A z < b. If z
% is not provided, the origin is assumed to be an interior point.
%
% IND_MINREP is a logical index array whose true elements correspond
% to the nonredundant halfplanes of the polytope {x | A x <= b}. The
% logical vector IND_NOT_VERIFIED refers to the indices of the
% inequalities that were not marked as nonredundant, but could not be
% verified as redundant.  For safety, these inequalities are also
% included in IND_MINREP (i.e., for difficult problems this function
% may return a nonminimal representation).
%
% IND_TO_CHECK can be used to only check a subset of the
% halfplanes. Its true elements represent the halfspaces (rows) to
% check for nonredundancy. This can be used as a containment check as
% follows:
%
% A12 = [A1; A2];
% b12 = [b1; b2];
% ind_to_check = false(size(b12));
% ind_to_check(1:numel(b1)) = true;
% [~, ~, ind_nred] = indicate_nonredundant_halfplanes([A1; A2], [b1; b2], ind_to_check);
% is_contained = not(any(ind_nred & ind_to_check));
%
% The rationale for the above calculations is that if
%
%   P1 = {x | A1 x <= b1}  and  P2 = {x | A2 x <= b2},
%
% and we want to test if P1 contains P2, we can test if all halfplanes
% of P1 are redundant with respect to P2.
%
% See also: active_set_solver, find_BFS
%
% % Examples
% % --------
%
% % The below ASCII plot depicts a 2d case where a polyhedron (marked by
% % *) is represented as a set of inequalities and x1<2 is a redundant
% % constraint.
% %
% %  ^ x2
% %  |      x1<1
% %  |        |   /|
% %  |        |  / |
% %  |        | /  |
% %  |        |/   |
% %  |        +    |
% %  |       /|    |
% %  x2<x1+1/*|    |
% %  |     /**|    |
% % 0+--- /***|    |<--x1<2 (redundant constraint)
% %  |   /****|    |
% % -+--+---+-+----+----- x2>-0.3
% %  | /    : |    |
% %  |/     : |    |
% %  +------+-+----+-------->
% % /       0 |    |         x1
%
% % This example is the same as the ASCII plot above
% Px1 = [-1  1;...
%         1  0;...
%         0 -1;...
%         1  0];
% Pc1 = [  1;...
%          2;...
%        0.3;...
%          1];
% 
% tic
% [Qx1, Qc1] = indicate_nonredundant_halfplanes(Px1, Pc1);
% solve_time1 = toc;
%
% % This is a somewhat more practically realistic example
% m = 5000;
% n = 10;
% Px2 = rand(m,n) - 0.5;
% Pc2 = rand(m,1) + 0.5;
%
% tic
% [Qx2, Qc2] = indicate_nonredundant_halfplanes(Px2, Pc2);
% solve_time2 = toc;
%
% % Print 
% fprintf('\n\n---------------------------------------------------------\n')
% fprintf('EXAMPLE\t\t\tRUNTIME\t\t\tELIMINATED INEQUALITIES\n')
% fprintf('---------------------------------------------------------\n')
% fprintf('Example 1\t\t%f sec\t\t%d\n',solve_time1, size(Pc1,1)-size(Qc1,1))
% fprintf('Example 2\t\t%f sec\t\t%d\n',solve_time2, size(Pc2,1)-size(Qc2,1))

if nargin < 3 || isempty(ind_to_check)
    ind_to_check = true(size(b));
end

%Check if the origin is an interior point
should_shift_back = false;
b_tol = 1e-10;
if any(b < b_tol)
    if nargin == 4 && any(z)
        %Shift set
        b = b - A*z;
        should_shift_back = true;
    else
        error('The origin is not an interior point or located too close to the boundary. Provide an interior point z.')
    end
    
    if any(b < b_tol)
        error('The input z is not an interior point or located too close to the boundary.')
    end
end


% Small modifications on halfplane description
[A, b] = normalize_halfplane_description(A, b);
[A, b, is_symmetric, Isym] = make_set_symmetric(A, b);

% Need to sort ind_to_check after having made set symmetric
ind_to_check = ind_to_check(Isym);
Iinvsym(Isym) = 1:length(Isym); % Quick formula for permutation inverse

% Remove duplicate halfplanes
[~,ind_notdup] = unique_tol([A,b], 1e-5); % Note: unique_tol, not uniquetol
n_notdup = find(ind_notdup);
ind_duplicate = ~ind_notdup;

% Find problem dimension
[m, n] = size(A);

% Calculate norms (which are used to warm-start the active set solver)
norms = rownorm(A);

%% Identify halfplanes

% Some logical arrays
ind_remain = ind_notdup & ind_to_check;    
ind_nred   = false(size(b));
ind_red    = false(size(b));
ind_failed = false(size(b));
ind_detected_nred    = false(size(b));

% Initialize
ind_active = false(size(b));
x = [];

MAX_IT = sum(ind_remain); % Maximum number of iterations
for it = 1:MAX_IT
    
    % Choose halfplane to identify
    if any(ind_remain)
        if it == 1 || isempty(x)
            j = find(ind_remain, 1);
        else
            % Calculate inner products
            inner_products = A * x;
            
            % Find angles between remaining halfplanes and solution x
            Ix_cos = (inner_products(ind_remain) ./ (norm(x) * norms(ind_remain)));
            [~, jj] = max(Ix_cos);
            find_ind_remain = find(ind_remain);
            j = find_ind_remain(jj);
        end
    else
        break
    end
    
    ind_notdup = ind_notdup & ~ind_red;  % reduce size of A, b sent to solver, to reduce cmoputational time for matrix-vector multiplication
    n_notdup = find(ind_notdup);
    
    [x,ind_active_notdup,ind_detected_nred_notdup,optimal_solution_found] = ...
         active_set_solver(-A(j,:)',A(ind_notdup, :),b(ind_notdup),x,ind_active(ind_notdup));

    % A and b are "squished" if duplicate halfplanes were found, so we
    % need to be careful when updating ind_active and ind_BFS;
    ind_active(:) = false;
    ind_active(n_notdup(ind_active_notdup)) = true;
    ind_detected_nred(:) = false;
    ind_detected_nred(n_notdup(ind_detected_nred_notdup)) = true;
    
    if ( ~ind_active(j) && optimal_solution_found && abs(A(j, :) * x - b(j)) >= 1e-12 )
        % Note: Tolerance for judging a halfplane as redundant should be conservatively set.
        ind_red(j) = true;
    else
        ind_nred(j) = true;
    end
    
    if ~optimal_solution_found
        ind_failed(j) = true;
    end
    
    % No matter what, j has been checked, so remove it
    ind_remain(j) = false;
    
    % Check if phase 1 has to be performed in next iteration
    if sum(ind_active) < n
        ind_active(:) = false;
        x = [];
    end
    
    ind_nred   = or(ind_nred, ind_detected_nred);
    ind_remain = and(ind_remain, ~ind_detected_nred);
    
    if is_symmetric
        % Halfplane mirrored to j
        ind_red(mod(j-1 + m/2, m) + 1) = ind_red(j);

        % Update nonredundant halfplanes
        ind_nred1 = ind_nred(1:m/2);
        ind_nred2 = ind_nred(m/2+1:end);
        ind_nred = or(ind_nred, [ind_nred2; ind_nred1]); % Note the reversed order

        ind_remain(mod(j-1 + m/2, m) + 1) = ind_remain(j);

        ind_failed(mod(j-1 + m/2, m) + 1) = ind_failed(j);
        
        % Update ind_remain
        ind_remain1 = ind_remain(1:m/2);
        ind_remain2 = ind_remain(m/2+1:end);
        ind_remain = and(ind_remain,[ind_remain2; ind_remain1]);
    end
    
end

% Return nonredundant and failed halfplanes
ind_minrep = or(ind_failed, ind_nred);
ind_not_verified = and(ind_failed,~ind_nred);

if is_symmetric
    assert(mod(sum(ind_minrep), 2) == 0, 'Uneven number of halfplanes returned');
    assert(mod(sum(ind_nred), 2) == 0, 'Uneven number of nonredundant halfplanes returned');
end

% Permute back index vectors (due to making set symmetrical)
ind_minrep = ind_minrep(Iinvsym);

% ind_minrep = ind_minrep(1:end-size(bbox, 1));
A = A(Iinvsym, :);
b = b(Iinvsym, :);
Amin = A(ind_minrep,:);
bmin = b(ind_minrep);

if should_shift_back
    bmin = bmin + Amin * z;
end


% end indicate_nonredundant_halfplanes
