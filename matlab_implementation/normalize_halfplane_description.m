function [H, h] = normalize_halfplane_description(H, h, normalize_H)
%NORMALIZE_HALFPLANE_DESCRIPTION  Normalize halfplane description
%
% [Hn, hn] = NORMALIZE_HALFPLANE_DESCRIPTION(H, h) divides each row in
% the halfplane description by the right-hand-side value if it is
% nonzero (to avoid zero-division) and non-inf (to avoid resulting
% NaNs).
%
% % Example
% % -------
%
% H = [ 1  1;...
%       0  1;...
%      -1  0;...
%       0 -1];
% h = [0;1;2;3];
%
% [Hn, hn] = normalize_halfplane_description(H, h);

if nargin < 3
    normalize_H = false;
end

if normalize_H

    rnH = rownorm(H);
    H = bsxfun(@rdivide, H, rnH);
    h = h./rnH;

else

    nz = ( (abs(h) > 1e-10) & not(isinf(h)) );
    hnz = h(nz);
    abs_hnz = abs(hnz); % avoid division with negative-signed rhs since it
                        % changes meaning of inequality

    H(nz, :) = bsxfun(@rdivide, H(nz,:), abs_hnz);
    h(nz, 1) = hnz./abs_hnz;

end

% end normalize_halfplane_description
