function B = learn_basis(X, S, Binit)
% Learning basis using Lagrange dual (with basis normalization)
%
% This code solves the following problem:
% 
%    minimize_B   0.5*||X - B*S||^2
%    subject to   ||B(:,j)||_2 <= l2norm, forall j=1...size(S,1)
% 
% The detail of the algorithm is described in the following paper:
% 'Efficient Sparse Coding Algorithms', Honglak Lee, Alexis Battle, Rajat Raina, Andrew Y. Ng, 
% Advances in Neural Information Processing Systems (NIPS) 19, 2007
%
% Written by Honglak Lee <hllee@cs.stanford.edu>
% Copyright 2007 by Honglak Lee, Alexis Battle, Rajat Raina, and Andrew Y. Ng

L = size(X,1);
N = size(X,2);
M = size(S, 1);

SSt = S*S';
XSt = X*S';

if exist('Binit', 'var')
    dual_lambda = diag(Binit\XSt - SSt);
else
    dual_lambda = 10*abs(rand(M,1)); % any arbitrary initialization should be ok.
end

c = 1;
trXXt = sum(sum(X.^2));

lb=zeros(size(dual_lambda));
options = optimset('GradObj','on', 'Hessian','on','Algorithm','trust-region-reflective');
%  options = optimset('GradObj','on', 'Hessian','on', 'TolFun', 1e-7);

% [x, fval, exitflag, output] = fmincon(@(x) fobj_basis_dual(x, SSt, XSt, X, c, trXXt), dual_lambda, [], [], [], [], lb, [], [], options);

for i = 1:10
    %newton's method
    [f,g,H] = fobj_basis_dual(dual_lambda,SSt,XSt,X,c,trXXt);
    %sing_val = svd(H);
    %reg = sing_val(floor(.7*length(sing_val)));
    reg = .5;
    delta_lambda = cgls(H,-g,dual_lambda,100,1e-6,reg^2);%pcg(H,-g);
    dual_lambda = dual_lambda + delta_lambda;
end

% dual_lambda= x;

% Bt = (SSt+diag(dual_lambda)) \ XSt';
% B_dual= Bt';

% output.iterations
%dual_lambda= x;

Bt = (SSt+diag(dual_lambda)) \ XSt';
B = Bt';

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,g,H] = fobj_basis_dual(dual_lambda, SSt, XSt, X, c, trXXt)
% Compute the objective function value at x
L= size(XSt,1);
M= length(dual_lambda);

SSt_inv = (SSt + diag(dual_lambda))\eye(size(SSt));%inv(SSt + diag(dual_lambda));

% trXXt = sum(sum(X.^2));
if L>M
    % (M*M)*((M*L)*(L*M)) => MLM + MMM = O(M^2(M+L))
    f = -trace(SSt_inv*(XSt'*XSt))+trXXt-c*sum(dual_lambda);
    
else
    % (L*M)*(M*M)*(M*L) => LMM + LML = O(LM(M+L))
    f = -trace(XSt*SSt_inv*XSt')+trXXt-c*sum(dual_lambda);
end
%f= -f;

if nargout > 1   % fun called with two output arguments
    % Gradient of the function evaluated at x
    g = zeros(M,1);
    temp = XSt*SSt_inv;
    g = sum(temp.^2) - c;
    g = g';
    %g= -g;
    
    
    if nargout > 2
        % Hessian evaluated at x
        % H = -2.*((SSt_inv*XSt'*XSt*SSt_inv).*SSt_inv);
        H = -2.*((temp'*temp).*SSt_inv);
        %H = -H;
    end
end

return