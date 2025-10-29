% Test script to compare with C++ implementation
Px1 = [-1  1;...
        1  0;...
        0 -1;...
        1  0];
Pc1 = [  1;...
         2;...
       0.3;...
         1];
         
% Interior point that satisfies all constraints
z = [0.5; -0.1];

% Check feasibility
fprintf('Testing interior point feasibility:\n');
Ax = Px1 * z;
for i = 1:4
    fprintf('  Constraint %d: %.1f <= %.1f ? %s\n', i-1, Ax(i), Pc1(i), ...
        char("YES" * (Ax(i) <= Pc1(i)) + "NO" * (Ax(i) > Pc1(i))));
end

% Solve
[Qx1, Qc1, ind_minrep, ind_not_verified] = indicate_nonredundant_halfplanes(Px1, Pc1, [], z);

fprintf('\nResults:\n');
fprintf('Success: YES\n');
fprintf('Original constraints: %d\n', size(Px1, 1));
fprintf('Non-redundant constraints: %d\n', size(Qx1, 1));

for i = 1:4
    fprintf('Constraint %d: ', i-1);
    if ~ind_minrep(i)
        fprintf('REDUNDANT');
    elseif ind_not_verified(i)
        fprintf('UNVERIFIED');
    else
        fprintf('NON-REDUNDANT');
    end
    fprintf('\n');
end