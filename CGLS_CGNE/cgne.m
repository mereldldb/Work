function [X, r, res1, res2, err, obj] = cgne(A, b, r0, P, C, R, m_iter, eps, lambda, x_mod, number)

[m,n] = size(A);

R_chol = chol(R);

r     = r0;

r_mod = C\(b-A*x_mod);

%r = zeros(length(r0),1);
if ( lambda > 0 )
   x  = (1/lambda)*(R\(A'*r));
else
   r = zeros(m,1);
   x = zeros(n,1);
end
s     = b-A*x-C*r;
z     = P\s;
p     = z;
Cp    = C*p;
q     = A'*p;
%Rq    = R\q;
Rq = R_chol\(R_chol'\q);

it = 1;
gamma = dot(s,z);
tol = eps*sqrt(gamma);

res1(1) = norm(r);
res2(1) = sqrt(gamma);
err(1)  = norm(x-x_mod);

X = [];
while it <= m_iter%( sqrt(gamma) > tol & it <= m_iter)
tic
   xi = dot(q,Rq) + lambda*dot(p,Cp);
   alpha = gamma/xi;
   r = r + (lambda*alpha)*p; 
   x = x + alpha*Rq; 
   s = s - alpha*(A*Rq+lambda*Cp);
   z = P\s;

   beta = 1/gamma;
   gamma = dot(s,z);
   beta = gamma*beta;
   p  = z + beta*p; Cp = C*p;
   q  = A'*p;
   %Rq = R\q;
   Rq = R_chol\(R_chol'\q);
   
   el_time(it) = toc;

   res1(it) = norm(r);
   %err2(it) = norm(r-r_mod);
   res2(it) = sqrt(gamma);
   err(it)  = norm(x-x_mod);
   obj(it) = 1/2*abs((A*x-b)'*(C\(A*x-b)))+1/2*lambda*abs(x'*R*x);

   if mod(it,number) == 0
       X = [X x];
   end   
   it = it+1;

end
avg_time = mean(el_time)