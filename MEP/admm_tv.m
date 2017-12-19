function [x,err] = admm_tv(A, b, x0, P, C, R, F, iter_admm, tol_admm, iter_CG, tol_CG, lambda, a, x_mod)
u = sparse(length(F*x0),1);
v = u;
it = 0;
dif = Inf;
it = 0;
err = norm( x0 - x_mod );
while (dif > tol_admm && it < iter_admm)
    it = it+1;
    x = cgls_admm(A, b, x0, P, C, R, F, iter_CG, tol_CG, a, x_mod, u, v);
    v = wthresh(F*x+u,'s',lambda/a);
    u = u+F*x-v;
    dif = norm(x - x0);
    x0 = x;
    err = [err, norm( x - x_mod )];
end
