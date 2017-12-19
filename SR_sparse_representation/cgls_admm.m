function [x, res1, res2, err] = cgls_admm(A, b, x0, m_iter, eps, a, u, v)

[m,n] = size(A);

F = eye(length(x0));
x     = x0; Rx = x;
r     = (b-A*x);
s     = A'*r - a*x + a*(v-u);
%z     = P\s;

%p     = z;
p = s;
Ap    = A*p;
CAp   = Ap;
Rp = p;

res1(1) = norm(r);
res2(1) = norm(s);
%err(1)  = norm(x-x_mod);

it = 1;
%gamma = dot(s,z);
gamma = dot(s,s);
tol = eps*sqrt(gamma);

while ( sqrt(gamma) > tol & it <= m_iter)

   xi = dot(Ap,CAp) + a*dot(p,Rp);
   alpha = gamma/xi;

   x = x + alpha*p; Rx = Rx + alpha*Rp;
   %r = r - alpha*CAp;
   r = (b-A*x);
   s = A'*r - a*Rx + a*(v-u);
   %z = P\s;

   beta = 1/gamma;
   %gamma = dot(s,z);
   gamma = dot(s,s);
   beta = gamma*beta;
   %p  = z + beta*p;
   p = s + beta*p;
   Ap = A*p;
   CAp = Ap;
   Rp = p;

   it = it+1;
   res1(it) = norm(r);
   res2(it) = norm(s);
   %err(it)  = norm(x-x_mod);

end
