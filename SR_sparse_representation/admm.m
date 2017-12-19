function [x,obj] = admm(A, b, x0, iter_admm, tol_admm, iter_CG, tol_CG, lambda, a)
u = sparse(length(x0),1);
v = u;
it = 0;
dif = Inf;
it = 0;
while (dif > tol_admm && it < iter_admm)
    it = it+1;
    x = cgls_admm(A, b, x0, iter_CG, tol_CG, a, u, v);
    v = wthresh(x+u,'s',lambda/a);
    u = u+x-v;
    dif = norm(x - x0);
    x0 = x;
    obj = norm(A*x - b)^2+lambda*norm(x,1);
end

end