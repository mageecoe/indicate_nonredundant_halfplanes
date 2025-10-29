function [x, ind_active] = find_BFS(A, b, z)
%FIND_BFS  Find a vertex of the set {x | Ax <= b}
%
% x = FIND_BFS(A, b) finds a vertex of the set {x | Ax <= b} by
% performing iterative rayshots starting at an interior point. To
% guarantee success with this syntax, the set needs to be compact and
% contain the origin.
%
% [x, IND_ACTIVE] = FIND_BFS(A, b, z) initializes x at the interior
% point z, and performes a rayshot in the direction s = -A(1,:)'. When
% an inequality is activated, x is updated and a new rayshot is
% performed in the null space of the active constraint. The procedure
% is repeated until a vertex has been found. IND_ACTIVE returns all
% constraints that are active at the vertex.
%
% Note that z has to be feasible, i.e. Az <= b has to be satisfied.
%
% See also: active_set_solver
%
% % Example
% % -------
%
% % Generate data
% m = 500;
% n = 5;
% A = rand(m,n) - 0.5;
% b = ones(m,1);
% 
% % Solve
% tic
% [x, ind_active] = find_BFS(A,b,zeros(n,1));
% solve_time_active_set =  toc;
% 
% % Print 
% fprintf('\nRun time: %f sec\n', solve_time_active_set)
%
%

% Find problem dimensions
[m,n] = size(A);

if nargin < 3 || not(any(z))
    if all(b > 0)
        z = zeros(n,1);
    else
        error(['The origin is not an interior point and no interior ' ...
               'point was provided. When the origin is not an ' ...
               'interior point an interior point, z, has to be provided.'])
    end

end

% Keep track of active set
ind_active = false(m,1);
n_active   = [];

% Initialize variables and search direction
s = -A(1,:)';
x = z;

% Find BFS
max_it = n;
for it=1:max_it

    % Find step size
    t = (b - A*x)./(A*s);
    t(ind_active)   = Inf;    % Avoid halfplanes already active
    t(t < 0)        = Inf;    % Get rid of solutions in negative direction while keeping size of t intact
    t(A*s <= 1e-12) = Inf;    % Make halfplanes that become "more feasible" along s permeable
    [min_t, n_cand] = min(t);

    if any(n_active == n_cand)
        error('Added already active halfplane. This indicates a problem with the algorithm.')
    end

    % Update active set
    ind_active(n_cand) = true;
    n_active(end+1) = n_cand;
    
    % Check stopping criterion
    if sum(ind_active) == n
        x = A(ind_active,:) \ b(ind_active);
        break
    end

    % Find QR-factorization
    [Q, R] = qr(A(ind_active,:)');
    m_active = sum(ind_active);
    R1 = R(1:m_active,:);
    Q1 = Q(:,1:m_active);
    Q2 = Q(:,m_active+1:end);
    
    if rcond(R(1:m_active,:)) < 1e-15
        warning('The problem is ill-conditioned. Consider using an alternative method to find an initial vertex.')
    end

    % Calculate x = arg min  0.5 || x - x_approx ||_2^2
    %                   s.t. A_active *x == b_active
    x_approx = x + min_t * s;
    res_approx = A(ind_active,:)*x_approx - b(ind_active);
    x = x_approx - Q1*( R1' \ res_approx );

    % Check feasibility
    if any(A(~ind_active,:)*x - b(~ind_active) > 0)
        error('Primal feasibility lost! This is probably due to an ill-conditioned A(ind_active,:).')
    end

    % Calculate new search direction
    s = Q2(:,end);

end

% end find_BFS
