function [x, res1, res2] = cgls(A, b, x0, m_iter, eps, lambda)
 
[m,n] = size(A);
 
%C_chol = chol(C);
 
x     = x0; Rx = x;
r     = (b-A*x);
%r     = C_chol\(C_chol'\(b-A*x));
s     = A'*r - lambda*x;
z     = s;
 
p     = z;
Ap    = A*p;
%CAp = C_chol\(C_chol'\Ap);
CAp   = Ap;
Rp = p;
 
res1(1) = norm(r);
res2(1) = norm(s);
 
it = 1;
gamma = dot(s,z);
tol = eps*sqrt(gamma);
 
while ( sqrt(gamma) > tol & it <= m_iter)
   tic
   xi = dot(Ap,CAp) + lambda*dot(p,Rp);
   alpha = gamma/xi;
 
   x = x + alpha*p; Rx = Rx + alpha*Rp;
   r = r - alpha*CAp;
   s = A'*r - lambda*Rx;
   z = s;
 
   beta = 1/gamma;
   gamma = dot(s,z);
   beta = gamma*beta;
   p  = z + beta*p;
   Ap = A*p;
   %CAp = C_chol\(C_chol'\Ap);
   CAp = Ap;
   Rp = p;
 
   it = it+1;
   res1(it) = norm(r);
   res2(it) = sqrt(gamma);
   el_time(it) = toc;    

    
end