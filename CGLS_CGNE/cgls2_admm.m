function [X, res1, res2, err] = cgls2_admm(A, b, x0, P, C, R, F, m_iter, eps, a, x_mod, u, v)

[m,n] = size(A);

x     = x0; Rx = R*x;
r     = C\(b-A*x);
s     = A'*r - a*R*x + a*F'*(v-u);
%z     = P\s;

%p     = z;
p = s;
Ap    = A*p;
CAp   = C\Ap;
Rp = R*p;

res1(1) = norm(r);
res2(1) = norm(s);
err(1)  = norm(x-x_mod);

it = 1;
%gamma = dot(s,z);
gamma = dot(s,s);
tol = eps*sqrt(gamma);
X = [];
while it <= m_iter%( sqrt(gamma) > tol & it <= m_iter)

   xi = dot(Ap,CAp) + a*dot(p,Rp);
   alpha = gamma/xi;

   x = x + alpha*p; Rx = Rx + alpha*Rp;
   %r = r - alpha*CAp;
   r = C\(b-A*x);
   s = A'*r - a*Rx + a*F'*(v-u);
   %z = P\s;

   beta = 1/gamma;
   %gamma = dot(s,z);
   gamma = dot(s,s);
   beta = gamma*beta;
   %p  = z + beta*p;
   p = s + beta*p;
   Ap = A*p;
   CAp = C\Ap;
   Rp = R*p;

   res1(it) = norm(r);
   res2(it) = norm(s);
   err(it)  = norm(x-x_mod);
   it = it+1;

   X = [X x];
end
