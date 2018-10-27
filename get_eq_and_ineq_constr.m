function [Hineq, hineq, Heq, heq] = get_eq_and_ineq_constr(H, h)
%GET_EQ_AND_INEQ_CONSTR  Partition constraints from inequalities to equalities and inequalities
%
% [Hineq, hineq, Heq, heq] = GET_EQ_AND_INEQ_CONSTR(H, h) identifies
% inequality pairs on the form H(i,:)*x<h(i), H(j,:)*x>h(j) and
% reduces each such pair to a single equality constraint,
% H(i,:)*x==h(i). The output is a partition of constraints into
% inequality and equality constraints.
%
% % Example:
% % -------
%
% n = 2;
% H = [eye(n);-eye(n)];
% h = [0;0;0;0];
% [Hineq, hineq, Heq, heq] = get_eq_and_ineq_constr(H, h);
%
% n = 3;
% H = [eye(n);-eye(n)];
% h = [0;2;3;0;1;1];
% [Hineq, hineq, Heq, heq] = get_eq_and_ineq_constr(H, h);

% Normalize (H,h)
[H, h] = normalize_halfplane_description(H, h);

n = size(H, 1);

%Hineq = zeros(0, size(H, 2));
%hineq = zeros(0, 1);
Hineq = zeros(n, size(H, 2));
hineq = zeros(n, 1);

%Heq = zeros(0, size(H, 2));
%heq = zeros(0, 1);
Heq = zeros(n, size(H, 2));
heq = zeros(n, 1);

% The following algorithm processes rows in H from top to bottom

k_eq = 0;
k_ineq = 0;
rows_left = true(size(h));
Hh = [H, h];
while any(rows_left);
    k = find(rows_left, 1);
    rows_left(k) = false;
    Hkhk = [H(k,:), h(k)];
    % bsxf = bsxfun(@plus, Hh, Hkhk);
    % iseq = sum(abs(bsxf),2)==0;
    iseq = all(bsxfun(@eq, Hh, -Hkhk), 2);
    if any(iseq)
        k_eq = k_eq + 1;
        Heq(k_eq, :) = H(k, :);
        heq(k_eq, 1) = h(k, 1);
        % Heq = [Heq; ...
        %        H(1, :)];
        % heq = [heq; ...
        %        h(1, :)];
        % % Remove rows that are duplicates
        % H(1 + find(iseq), :) = [];
        % h(1 + find(iseq), :) = [];
        rows_left(iseq) = false;
    else
        % Add row to ineq part
        k_ineq = k_ineq + 1;
        Hineq(k_ineq, :) = H(k, :);
        hineq(k_ineq, 1) = h(k, 1);
        % Hineq = [Hineq; ...
        %          H(1, :)];
        % hineq = [hineq; ...
        %          h(1)];
    end
    % % Remove processed row
    % H(1, :) = [];
    % h(1, :) = [];
end

% % Add last unprocessed row (if it was not a duplicate row)
% Hineq = [Hineq; ...
%          H];
% hineq = [hineq; ...
%          h];

Heq(k_eq+1:end, :) = [];
heq(k_eq+1:end, :) = [];
Hineq(k_ineq+1:end, :) = [];
hineq(k_ineq+1:end, :) = [];


% Fix dimensions if empty arrays
if ( isempty(Hineq) && isempty(hineq) )
    Hineq = zeros(0, size(H, 2));
    hineq = zeros(0, size(h, 2));
end

% end get_eq_and_ineq_constr
