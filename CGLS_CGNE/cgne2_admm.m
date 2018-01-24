function [X, r, res1, res2, err] = cgne2_admm(A, b, r0, P, C, R, F, m_iter, eps, lambda, x_mod, u, v)

[m,n] = size(A);

R_chol = chol(R);

r     = r0;

r_mod = C\(b-A*x_mod);

%r = zeros(length(r0),1);
if ( lambda > 0 )
   y = 1/lambda*r;
   x  = 1/lambda*(R\((A'*r)+lambda*(F'*(v-u))));
else
   r = zeros(m,1);
   x = zeros(n,1);
   y = zeros(m,1);
end
s     = b-C*r-A*x;
%z     = P\s;
p     = s;
Cp    = C*p;
q     = A'*p;
Rq    = R\q;
%Rq = R_chol\(R_chol'\q);

it = 1;
gamma = dot(s,s);
tol = eps*sqrt(gamma);

X = [];
while it <= m_iter %( sqrt(gamma) > tol & it <= m_iter )
   tic 
   xi = dot(q,Rq) + lambda*dot(p,Cp);
%    xi = 1/lambda*dot(q,Rq) + dot(p,Cp);
   alpha = gamma/xi;
   y = y + alpha*p;
%    r = r + alpha*p;
   x = R\((A'*y)+(F'*(v-u)));

%    x = 1/lambda*R\((A'*r)+lambda*(F'*(v-u)));
%   x = x + alpha*Rq;% + R\(F'*(u-v));
%    s = b-A*x-C*r;
%    s = s - alpha*(Cp+1/lambda*A*Rq);
    s = s - alpha*(lambda*Cp+A*Rq);
%   z = P\s;

   beta = 1/gamma;
   gamma = dot(s,s);
   beta = gamma*beta;
   p  = s + beta*p; Cp = C*p;
   q  = A'*p;
   Rq = R\q;
   %Rq = R_chol\(R_chol'\q);
   
   el_time(it) = toc;

   res1(it) = norm(r);
   %err2(it) = norm(r-r_mod);
   res2(it) = sqrt(gamma);
   err(it)  = norm(x-x_mod);
   %obj(it) = real(norm(A*x-b)^2+1/2*lambda*x'*R*x);

   it = it+1;

   X = [X x];
end
%avg_time = mean(el_time)