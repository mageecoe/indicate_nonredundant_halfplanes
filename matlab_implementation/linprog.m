## linprog.m  — simple drop-in replacement using glpk()

function [x, fval, exitflag, output] = linprog(f, A, b, Aeq, beq, lb, ub)
  % LINPROG  Linear programming via GLPK
  %
  % Syntax (compatible subset):
  %   [x, fval] = linprog(f, A, b)
  %   [x, fval] = linprog(f, A, b, Aeq, beq, lb, ub)
  %
  % Note: this function uses Octave's built-in GLPK solver.

  if nargin < 7, ub = []; end
  if nargin < 6, lb = []; end
  if nargin < 5, beq = []; end
  if nargin < 4, Aeq = []; end

  % Combine inequality and equality constraints
  A_all = A;
  b_all = b;
  ctype = repmat("U", size(A,1), 1);  % U = upper bound (A*x <= b)

  if ~isempty(Aeq)
    A_all = [A_all; Aeq];
    b_all = [b_all; beq];
    ctype = [ctype; repmat("S", size(Aeq,1), 1)]; % S = equality (Aeq*x = beq)
  end

  % Default lower/upper bounds
  n = numel(f);
  if isempty(lb), lb = -inf(n,1); end
  if isempty(ub), ub =  inf(n,1); end

  % Run solver
  sense = 1; % 1 = minimize
  [x, fval, status, extra] = glpk(f, A_all, b_all, lb, ub, ctype, [], sense);

  % Convert GLPK status to linprog-style flag
  exitflag = 1; % assume success
  if status ~= 0
    exitflag = 0;
  end

  output.algorithm = "glpk";
  output.status = status;
  output.time = extra.time;
end
