function [H, h, B] = subspacify_polyhedron(H, h);
%SUBSPACIFY_POLYHEDRON  Represent polyhedron in subspace
%
% [Hout, hout, B] = SUBSPACIFY_POLYHEDRON(Hin, hin) calculates a
% representation (Hout,hout) of the polyhedron (Hin, hin) projected
% onto the subspace that is spanned by the principal axes of the
% polyhedron (Hin, hin). A basis B is returned for the subspace.
%
% See also: superspacify_polyhedron
%
% % Example
% % -------
%
% H = [eye(2); -eye(2)];
% h = [1; 0; 1; 0];
%
% [Hs, hs, B] = subspacify_polyhedron(H, h);

[Hineq, hineq, Heq, heq] = get_eq_and_ineq_constr(H, h);
B = null(Heq);
H = (B \ Hineq')';
h = hineq;

% end subspacify_polyhedron
