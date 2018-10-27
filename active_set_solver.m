function [x, ind_active, ind_detected_nred, optimal_solution_found] = active_set_solver(f, A, b, x, ind_active)
%ACTIVE_SET_SOLVER  Linear program solver
%
% x = ACTIVE_SET_SOLVER(f, A, b) solves the linear program
%
%  min  f'x
%  s.t. Ax \leq b
%
% using a primal active set method. The method is tailored for use in
% the external function indicate_nonredundant_halfplanes. For more
% details about the algorithm, see the paper
%
% Emil Klintberg, Magnus Nilsson, Lars Johannesson Mårdh. "A primal
% active-set minimal representation algorithm for polytopes with
% application to invariant-set calculations".
%
% The general syntax of the function is
%
% [x, IND_ACTIVE, ...
%     IND_DETECTED_NRED, ...
%     OPTIMAL_SOLUTION_FOUND] = ACTIVE_SET_SOLVER(f, A, b, x, IND_ACTIVE)
%
% The active-set method benefits from warm-starting, so an inital
% active set can be specified using IND_ACTIVE. The optimal active set
% is in turn returned in IND_ACTIVE. If an initial active set is not
% provided, a starting point is searched for with find_BFS.
%
% The input x should be the point that is defined by the input
% IND_ACTIVE. Similarly, the output x is the point that is defined by
% the output IND_ACTIVE. This is always the case when this function is
% called by indicate_nonredundant_halfplanes.
%
% IND_DETECTED_NRED returns all constraints that were part of any
% active set over the iterations. This is exploited in
% indicate_nonredundant_halfplanes.
%
% OPTIMAL_SOLUTION_FOUND returns true if the optimal solution was
% found and false otherwise.
%
%
%
% See also: find_BFS
%
% % Example
% % -------
%
% % Generate data
% m = 2000;
% n = 10;
% f = rand(n,1) - 0.5;
% A = rand(m,n) - 0.5;
% b = ones(m,1);
% 
% % Solve
% tic
% [x_active_set, ...
%  ind_active, ...
%  ind_detected_nred, ...
%  optimal_solution_found] = active_set_solver(f, A, b);
% solve_time_active_set =  toc;
% assert(optimal_solution_found, 'Optimal solution not found!')
% 
% tic
% x_linprog = linprog(f, A, b);
% solve_time_linprog = toc;
% 
% % Print 
% fprintf('\n\n-----------------------------------------------------\n')
% fprintf('METHOD\t\t\tRUNTIME\t\t\tOPTIMAL VALUE\n')
% fprintf('-----------------------------------------------------\n')
% fprintf('Active set\t\t%f sec\t\t%f\n',solve_time_active_set, f'*x_active_set)
% fprintf('linprog\t\t\t%f sec\t\t%f\n',solve_time_linprog, f'*x_linprog)


% Find dimensions of problem
[~, n] = size(A);


% Check if phase 1 has to be performed
if nargin < 5 || not(any(ind_active))
    [x, ind_active] = find_BFS(A,b);
end

%Initialize
ind_detected_nred      = false(size(b));
optimal_solution_found = false;

n_active = find(ind_active);
I = eye(n);
I_reduced = eye(n-1);

%% Solve
MAX_IT = 10000;
for it = 1:MAX_IT

    % Add ind_active to ind_detected_nred (used in indicate_nonredundant_halfplanes)
    ind_detected_nred = or(ind_detected_nred,ind_active);

    % Find QR-factorization in first iteration
    if it == 1
        A_active = A(ind_active,:);
        b_active = b(ind_active);
        [Q,R] = qr(A_active');
    end

    % Calculate dual variables using QR factors
    Qf = -Q'*f;
    mu = R\Qf;
    %mu = -R \ (Q'*f);

    % Check stopping criterion
    if all(mu >= -1e-6)
        optimal_solution_found = true;
        break
    end

    %----------------------------------------------------------------------
    % Find and remove halfplane from active set
    %----------------------------------------------------------------------
    [~,n_remove] = min(mu);
    ind_active(n_active(n_remove)) = false;

    a_remove = A(n_active(n_remove),:);

    %----------------------------------------------------------------------
    % Project (negated) normal vector of removed halfplane onto nullspace
    %----------------------------------------------------------------------
    if it==1
        % Find A_reduced
        A_reduced = A(ind_active,:);

        % Find indices of halfplanes in A_reduced
        n_reduced = find(ind_active);

        % Find QR-factorization
        [Q_reduced,R_reduced] = qr(A_reduced');
    else
        % If A_reduced has changed
        if any(n_reduced == n_active(n_remove))
            % Find the index that has changed
            ind = find(n_reduced == n_active(n_remove));

            % Update A_reduced
            a_reduced_remove = A_reduced(ind,:);
            A_reduced(ind,:) = a_add;

            % Update vector of indices
            n_reduced(ind) = n_add;

            % Update QR-factorization
            [Q_reduced,R_reduced] = qrupdate(Q_reduced,R_reduced,a_add'-a_reduced_remove',I_reduced(:,ind));

        end
    end

    % Calculate projection
    s = Q_reduced(:, end);

    % Check that we are moving in the right direction
    if a_remove * s > 0
        s = -s;
    end

    %----------------------------------------------------------------------
    % Find maximum step size to remain primal feasible
    %----------------------------------------------------------------------

    % Find step sizes to activate inactive halfplanes
    p1 = b - A*x;
    p2 = A*s;

    step_sizes = p1./p2;
    step_sizes(p2 <= 1e-12 | ind_active) = Inf;

    %----------------------------------------------------------------------
    % Find variable to add and step size
    %----------------------------------------------------------------------
    [t,n_add] = min(step_sizes);

    % Perform some checks
    if isinf(t)
        error('The step size has infinite length. Assert that your provided set is compact and that your starting point z (the origin if not provided as input) is an interior point of the set.')
    end

    if any(n_add == find(ind_active))
        error('Trying to add an already active constraint. This indicates a problem with the algorithm.')
    end

    %----------------------------------------------------------------------
    % Evaluate step and update variables
    %----------------------------------------------------------------------

    % Update active constraints
    a_add = A(n_add,:);
    b_add = b(n_add);
    A_active(n_remove,:) = a_add;
    b_active(n_remove,:) = b_add;
    n_active(n_remove) = n_add;
    ind_active(n_add)  = true;

    % Update QR-factorization
    [Q, R] = qrupdate(Q, R, a_add'-a_remove', I(:,n_remove));

    %----- Update variables -----
    % Quit if matrix is close to singular
    if rcond(R) < 1e-15
        optimal_solution_found = false;
        ind_active = false(size(b));
        return
    end

    % Update variables
    x = Q * (R' \ b_active);

    % Correct variables if inaccurate
    if max(abs(A_active*x - b_active)) > 1e-9

        % Quit if matrix is close to singular
        if rcond(A_active) < 1e-15
            optimal_solution_found = false;
            ind_active = false(size(b));
            return
        end

        x = A_active\b_active;
    end

    % Verify that variables are feasible
    max_res = max(abs(A_active*x - b_active));
    if max(A*x-b) > max_res + 1e-7
        optimal_solution_found = false;
        ind_active = false(size(b));
        return
    end

    %----------------------------------------------------------------------
    % Check if maximum number of iterations has been reached
    %----------------------------------------------------------------------
    if it==MAX_IT
        warning('Maximum number of iterations reached!\n')
        optimal_solution_found = false;
        ind_active = false(size(b));
        return
    end

end

% end active_set_solver
