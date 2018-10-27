function [H, h, varargout] = superspacify_polyhedron(H, h, B);
%SUPERSPACIFY_POLYHEDRON  Represent polyhedron in superspace
%
% [Hout, hout] = SUPERSPACIFY_POLYHEDRON(Hin, hin, B) calculates a
% halfplane representation (Hout,hout) of the polyhedron (Hin, hin)
% lifted onto the superspace that the basis B lies in.
%
% See also: subspacify_polyhedron
%
% % Example
% % -------
%
% H = [-1;1];
% h = [1;1];
% B = [-1;0];
%
% [Hs, hs] = superspacify_polyhedron(H, h, B);

H * B';
H = H * B';
nullB = null(B');
m = size(nullB, 2);
if nargout <= 2
    % Represent equality constraints in H, h
    H(end+[1:2*m], :) = [nullB';-nullB'];
    h(end+[1:2*m], 1) = 0;
elseif nargout == 4
    Heq([1:2*m], :) = [nullB';-nullB'];
    heq([1:2*m], 1) = 0;
    varargout{1} = Heq;
    varargout{2} = heq;
else
    error('Bad number of output arguments in superspacify_polyhedron')
end


% end superspacify_polyhedron
