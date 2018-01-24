function [x,enrm,obj2] = admm_tv(A, b, x0, P, C, R, F, iter_admm, tol_admm, iter_CG, tol_CG, lambda, a, x_mod)
u = sparse(length(F*x0),1);
v = u;
dif = Inf;
it = 0;
enrm = [];
obj2 = [];
while it < iter_admm%(dif > tol_admm && it < iter_admm)
    it = it+1;
    [X, res1, ~, err] = cgls2_admm(A, b, x0, P, C, R, F, iter_CG, tol_CG, a, x_mod, u, v);
    x = X(:,end);
    v = wthresh(F*x+u,'s',lambda/a);
    u = u+F*x-v;
    dif = norm(x - x0);
    x0 = x;
    enrm = [enrm err];
    for k = 1:size(X,2)
        obj(k) = 1/2*abs((b-A*X(:,k))'*(C\(b-A*X(:,k))))+1/2*lambda*sum(abs(F*X(:,k)));
    end
    obj2 = [obj2 obj];
end
