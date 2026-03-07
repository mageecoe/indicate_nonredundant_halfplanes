% Comprehensive test comparing MATLAB, C++, and qhull implementations
fprintf('=== Polytope Redundancy Removal: Implementation Comparison ===\n\n');

%% Test Case 1: Simple 2D example
fprintf('Test Case 1: Simple 2D polytope\n');
fprintf('--------------------------------\n');

Px1 = [-1  1;...
        1  0;...
        0 -1;...
        1  0];
Pc1 = [  1;...
         2;...
       0.3;...
         1];

% Interior point
z = [0.0; 0.0];

% Test all implementations
fprintf('Running MATLAB implementation...\n');
[Qx1_matlab, Qc1_matlab, ind_minrep_matlab, ind_not_verified_matlab] = ...
    indicate_nonredundant_halfplanes(Px1, Pc1, [], z);

fprintf('Running qhull implementation...\n');
try
    [Qx1_qhull, Qc1_qhull, ind_minrep_qhull, ind_not_verified_qhull] = ...
        qhull_minimal_representation(Px1, Pc1, z);
    qhull_available = true;
catch ME
    fprintf('Qhull failed: %s\n', ME.message);
    qhull_available = false;
end

% Compare results
fprintf('\nResults:\n');
fprintf('Original constraints: %d\n', size(Px1, 1));
fprintf('MATLAB non-redundant: %d\n', size(Qx1_matlab, 1));
if qhull_available
    fprintf('Qhull non-redundant:  %d\n', size(Qx1_qhull, 1));
end

fprintf('\nConstraint-by-constraint comparison:\n');
for i = 1:4
    fprintf('Constraint %d: ', i-1);
    if ~ind_minrep_matlab(i)
        fprintf('MATLAB=REDUNDANT ');
    elseif ind_not_verified_matlab(i)
        fprintf('MATLAB=UNVERIFIED ');
    else
        fprintf('MATLAB=NON-REDUNDANT ');
    end
    
    if qhull_available
        if ~ind_minrep_qhull(i)
            fprintf('QHULL=REDUNDANT');
        else
            fprintf('QHULL=NON-REDUNDANT');
        end
    end
    fprintf('\n');
end

%% Test Case 2: 3D cube with redundant constraints
fprintf('\n\nTest Case 2: 3D cube with redundant constraints\n');
fprintf('------------------------------------------------\n');

% Standard cube: -1 <= x,y,z <= 1 plus redundant constraints
A_cube = [ 1  0  0;   % x <= 1
          -1  0  0;   % -x <= 1  
           0  1  0;   % y <= 1
           0 -1  0;   % -y <= 1
           0  0  1;   % z <= 1
           0  0 -1;   % -z <= 1
           1  0  0;   % x <= 2 (redundant)
           0  1  0];  % y <= 1.5 (redundant)

b_cube = [1; 1; 1; 1; 1; 1; 2; 1.5];
z_cube = [0; 0; 0];

fprintf('Running MATLAB implementation...\n');
[A_min_matlab, b_min_matlab, ind_minrep_matlab, ~] = ...
    indicate_nonredundant_halfplanes(A_cube, b_cube, [], z_cube);

if qhull_available
    fprintf('Running qhull implementation...\n');
    try
        [A_min_qhull, b_min_qhull, ind_minrep_qhull, ~] = ...
            qhull_minimal_representation(A_cube, b_cube, z_cube);
    catch ME
        fprintf('Qhull failed: %s\n', ME.message);
        qhull_available = false;
    end
end

fprintf('\nResults:\n');
fprintf('Original constraints: %d\n', size(A_cube, 1));
fprintf('MATLAB non-redundant: %d\n', size(A_min_matlab, 1));
if qhull_available
    fprintf('Qhull non-redundant:  %d\n', size(A_min_qhull, 1));
end

fprintf('\nExpected: Constraints 6 and 7 should be redundant\n');
fprintf('MATLAB: Constraint 6 %s, Constraint 7 %s\n', ...
    char("redundant" * ~ind_minrep_matlab(7) + "non-redundant" * ind_minrep_matlab(7)), ...
    char("redundant" * ~ind_minrep_matlab(8) + "non-redundant" * ind_minrep_matlab(8)));

if qhull_available
    fprintf('Qhull:  Constraint 6 %s, Constraint 7 %s\n', ...
        char("redundant" * ~ind_minrep_qhull(7) + "non-redundant" * ind_minrep_qhull(7)), ...
        char("redundant" * ~ind_minrep_qhull(8) + "non-redundant" * ind_minrep_qhull(8)));
end

%% Test Case 3: Random polytope
fprintf('\n\nTest Case 3: Random polytope\n');
fprintf('-----------------------------\n');

rng(42); % For reproducibility
m = 50;  % constraints
n = 5;   % dimension

A_rand = randn(m, n);
b_rand = ones(m, 1) + 0.5 * rand(m, 1);  % Ensure polytope is feasible
z_rand = zeros(n, 1);

% Check if origin is interior
origin_feasible = all(A_rand * z_rand < b_rand - 1e-6);
if ~origin_feasible
    z_rand = 0.1 * ones(n, 1);
end

fprintf('Running MATLAB implementation...\n');
tic;
[A_min_matlab, b_min_matlab, ind_minrep_matlab, ~] = ...
    indicate_nonredundant_halfplanes(A_rand, b_rand, [], z_rand);
matlab_time = toc;

if qhull_available
    fprintf('Running qhull implementation...\n');
    tic;
    try
        [A_min_qhull, b_min_qhull, ind_minrep_qhull, ~] = ...
            qhull_minimal_representation(A_rand, b_rand, z_rand);
        qhull_time = toc;
    catch ME
        fprintf('Qhull failed: %s\n', ME.message);
        qhull_available = false;
        qhull_time = NaN;
    end
end

fprintf('\nResults:\n');
fprintf('Original constraints: %d\n', m);
fprintf('MATLAB non-redundant: %d (%.2f ms)\n', size(A_min_matlab, 1), matlab_time * 1000);
if qhull_available
    fprintf('Qhull non-redundant:  %d (%.2f ms)\n', size(A_min_qhull, 1), qhull_time * 1000);
    
    % Check agreement
    if size(A_min_matlab, 1) == size(A_min_qhull, 1)
        fprintf('✓ Both methods found the same number of non-redundant constraints\n');
    else
        fprintf('⚠ Methods disagree on number of non-redundant constraints\n');
        fprintf('  Difference: %d constraints\n', abs(size(A_min_matlab, 1) - size(A_min_qhull, 1)));
    end
end

%% Summary
fprintf('\n\n=== SUMMARY ===\n');
if qhull_available
    fprintf('All implementations completed successfully.\n');
    fprintf('Qhull provides an independent reference for validation.\n');
else
    fprintf('Qhull not available - only MATLAB implementation tested.\n');
    fprintf('Install qhull for comprehensive comparison.\n');
end

fprintf('\nNote: Small differences between methods may occur due to:\n');
fprintf('- Different numerical tolerances\n');
fprintf('- Different tie-breaking in degenerate cases\n');
fprintf('- Numerical precision differences\n');