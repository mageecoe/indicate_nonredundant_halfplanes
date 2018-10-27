function out = rownorm(A, n)
%ROWNORM  Norm of columns of matrix
%
% N = ROWNORM(A) returns a column vector where each element represents
% the vector norm of each row of A.
%
% % Example
% % -------
%
% A = [3 1;...
%      4 0];
% b = rownorm(A)
%
% % Timing test
%
% A = rand(50000,6);
% tstart = tic;
% rownorm(A);
% tstop = toc(tstart)

if nargin < 2
    n = 2; % 2-norm by default
end

sizeA1 = size(A, 1);
out = zeros(sizeA1, 1);

%out = sum(A.^n, 2).^(1/n);

for idx = 1:sizeA1
    out(idx, 1) = norm(A(idx, :), n);
%     % out(idx, 1) = sqrt(A(idx,:)*A(idx,:)');
end

% end rownorm
