function [Amin, bmin, ind_minrep, ind_not_verified] = qhull_minimal_representation(A, b, interior_point)
%QHULL_MINIMAL_REPRESENTATION  Qhull-based minimal H-representation
%
% [Amin, bmin] = qhull_minimal_representation(A, b) returns a minimal
% representation of the polytope {x | A x <= b} using qhull's qhalf program.
%
% [Amin, bmin, IND_MINREP, IND_NOT_VERIFIED] = qhull_minimal_representation(A, b, interior_point)
%
% This function serves as a reference implementation and comparison tool
% for the custom active-set algorithm. It requires qhull to be installed
% and qhalf to be available in the system PATH.
%
% Inputs:
%   A - Constraint matrix (m x n)
%   b - Constraint vector (m x 1)  
%   interior_point - Optional interior point (n x 1). If not provided,
%                   the origin is assumed.
%
% Outputs:
%   Amin - Minimal constraint matrix
%   bmin - Minimal constraint vector
%   ind_minrep - Logical vector indicating non-redundant constraints
%   ind_not_verified - Always empty (qhull doesn't have unverified constraints)
%
% See also: indicate_nonredundant_halfplanes

if nargin < 3 || isempty(interior_point)
    interior_point = zeros(size(A, 2), 1);
end

[m, n] = size(A);

% Initialize outputs
ind_minrep = false(m, 1);
ind_not_verified = false(m, 1);

if m == 0
    Amin = A;
    bmin = b;
    return;
end

% Check if qhalf is available
[status, ~] = system('which qhalf > /dev/null 2>&1');
if status ~= 0
    error('qhalf not found in PATH. Please install qhull.');
end

% Validate interior point
for i = 1:m
    if A(i, :) * interior_point >= b(i) - 1e-10
        error('Interior point is not in the interior of the polytope');
    end
end

try
    % Create temporary input file for qhalf
    temp_file = tempname;
    fid = fopen(temp_file, 'w');
    
    if fid == -1
        error('Failed to create temporary file');
    end
    
    % Write qhalf input format
    % First line: dimension 1 interior_point
    fprintf(fid, '%d 1', n);
    for j = 1:n
        fprintf(fid, ' %.16g', interior_point(j));
    end
    fprintf(fid, '\n');
    
    % Second line: (dimension+1) number_of_halfspaces  
    fprintf(fid, '%d %d\n', n+1, m);
    
    % Halfspace definitions: coefficients offset
    % qhull uses ax + by + c <= 0, we have Ax <= b
    % Convert: Ax <= b becomes -Ax + b <= 0
    for i = 1:m
        for j = 1:n
            fprintf(fid, '%.16g ', -A(i, j));
        end
        fprintf(fid, '%.16g\n', b(i));
    end
    
    fclose(fid);
    
    % Run qhalf with Fx option to get non-redundant halfspace indices
    command = sprintf('qhalf Fx < %s', temp_file);
    [status, output] = system(command);
    
    % Clean up temporary file
    delete(temp_file);
    
    if status ~= 0
        error('qhalf failed with status %d', status);
    end
    
    % Parse output to get non-redundant indices
    lines = strsplit(strtrim(output), '\n');
    nonredundant_indices = [];
    
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if ~isempty(line) && ~isnan(str2double(line))
            idx = str2double(line);
            if idx >= 0 && idx < m
                nonredundant_indices = [nonredundant_indices, idx + 1]; % Convert to 1-based indexing
            end
        end
    end
    
    % Mark non-redundant constraints
    ind_minrep(nonredundant_indices) = true;
    
    % Build minimal representation
    Amin = A(ind_minrep, :);
    bmin = b(ind_minrep);
    
catch ME
    warning('qhull_minimal_representation failed: %s', ME.message);
    % Return original polytope on failure
    Amin = A;
    bmin = b;
    ind_minrep = true(m, 1);
    ind_not_verified = false(m, 1);
end

end