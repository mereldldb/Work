function [x,enrm,obj2] = admm_tv_2(A, b, r0, P, C, R, F, iter_admm, tol_admm, iter_CG, tol_CG, lambda, a, x_mod)
u = sparse(size(F,1),1);
v = u;
dif = Inf;
it = 0;
enrm = [];
obj2 = [];
x0 = zeros(size(A,2),1);
while (dif > tol_admm && it < iter_admm)
    it = it+1;
    [X, r0, ~, ~, err] = cgne2_admm(A, b, r0, P, C, R, F, iter_CG, tol_CG, a, x_mod, u, v);
%     y = pcg(a*C+A*(R\(A')),b-A*(R\(F'*(u-v))));
%     x = R\((A'*y)+(F'*(v-u)));
    x = X(:,end);
    v = wthresh(F*x+u,'s',lambda/a);
    u = u+F*x-v;
    %dif = norm(x - x0);
    x0 = x;
    enrm = [enrm err];
    for k = 1:size(X,2)
        obj(k) = 1/2*norm(b-A*X(:,k))^2+1/2*lambda*sum(abs(F*X(:,k)));
    end
    obj2 = [obj2 obj];
end