function [H, h, is_symmetric, I] = make_set_symmetric(H, h, tol)
%MAKE_SET_SYMMETRIC  Reorder rows according to convention if set is symmetric
%
% [H, h, is_symmetric] = MAKE_SET_SYMMETRIC(H, h) reorders rows in {H
% x <= h} according to convention if the set is symmetric.
%
% The reordering convention is to divide the representation into an
% upper and lower part, where the lower part contains the mirrored
% halfplanes
%
% % Example
% % -------
%
% % Create a random symmetric polytope
% m = 10000;
% n = 20;
% Hu = rand(m, n) - 0.5;
% H = [Hu; -Hu];
% h = ones(size(H, 1), 1);
% rowperm = randperm(size(H, 1));
% Hp = H(rowperm, :); hp = h(rowperm);
% tic;
% [Hs, hs, is_symmetric] = make_set_symmetric(Hp, hp);
% toc

if nargin < 3
    % numerical tolerance for equality
    tol = 1e-6;
end

is_symmetric = false;
nrows = numel(h);
I = [1:nrows]';
if mod(nrows, 2) == 0  % even number of rows
    [H, h] = normalize_halfplane_description(H, h);
    Hh = [H, h];
    % 1. sort rows
    [Hhs, J] = sortrows(Hh);
    % 2. flip lower part
    upper = [1:nrows/2];
    lower = [nrows/2+1:nrows];
    Hhs(lower,:) = flipud(Hhs(lower, :));
    J(lower) = flipud(J(lower));
    % 3. Check if lower part is symmetric with upper part
    Hs_upper = Hhs(upper, 1:end-1);
    hs_upper = Hhs(upper, end);
    Hs_lower = Hhs(lower, 1:end-1);
    hs_lower = Hhs(lower, end);
    is_symmetric = ( all(all(abs([Hs_upper, hs_upper] - [-Hs_lower, hs_lower]) < tol)) );
    if is_symmetric
        % Force set to be numerically symmetric
        I = J; % Only send permutation information if H, h is actually permuted.
        H = [Hs_upper; -Hs_upper];
        h = [hs_upper; hs_upper];
    end
end

% end make_set_symmetric
