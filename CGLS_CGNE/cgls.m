function [X, res1, res2, err, obj] = cgls(A, b, x0, P, C, R, m_iter, eps, lambda, x_mod, number)

[m,n] = size(A);

%C_chol = chol(C);

x     = x0; Rx = R*x;
r     = C\(b-A*x);
%r     = C_chol\(C_chol'\(b-A*x));
s     = A'*r - lambda*Rx;
z     = P\s;

p     = z;
Ap    = A*p;
%CAp = C_chol\(C_chol'\Ap);
CAp   = C\Ap;
Rp = R*p;

res1(1) = norm(r);
res2(1) = norm(s);
err(1)  = norm(x-x_mod);

it = 1;
gamma = dot(s,z);
tol = eps*sqrt(gamma);

X = [];
while it <= m_iter%( sqrt(gamma) > tol & it <= m_iter)
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
   %CAp = C_chol\(C_chol'\Ap);
   CAp = C\Ap;
   Rp = R*p;
   
   el_time(it) = toc;

   res1(it) = norm(r);
   res2(it) = sqrt(gamma);
   err(it)  = norm(x-x_mod);
   obj(it) = 1/2*abs((A*x-b)'*(C\(A*x-b)))+1/2*abs(lambda*x'*R*x);

   if mod(it,number) == 0
       X = [X x];
   end
   it = it+1;

end
avg_time = mean(el_time)