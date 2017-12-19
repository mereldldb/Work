function [x, res1, res2, err] = cgls2(A, b, x0, P, C, R, m_iter, eps, lambda, x_mod)

[m,n] = size(A);

x     = x0; Rx = R*x;
r     = C\(b-A*x);
s     = A'*r - lambda*Rx;
z     = P\s;

p     = z;
Ap    = A*p;
CAp   = C\Ap;
Rp = R*p;

res1(1) = norm(r);
res2(1) = norm(s);
err(1)  = norm(x-x_mod);

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
   z = P\s;

   beta = 1/gamma;
   gamma = dot(s,z);
   beta = gamma*beta;
   p  = z + beta*p;
   Ap = A*p;
   CAp = C\Ap;
   Rp = R*p;

   it = it+1;
   res1(it) = norm(r);
   res2(it) = sqrt(gamma);
   err(it)  = norm(x-x_mod);
   toc
end
