function [S,S2] = learn_sparsity(X,B,S_init,lambda)
    tic
    for i = 1:length(S_init(1,:))
        s_init = S_init(:,i);
        x = X(:,i);
%         tic
        [s,FitInfo] = lassoglm(B,x,'normal','Lambda',lambda);
%         disp('LASSO')
%         toc
%         tic
%         [s,obj] = admm(B,x,s_init,40,1e-10,10,1e-10,lambda,10*lambda);
        S(:,i) = s;
%         disp('ADMM')
%         toc
%         tic
%         disp('Feature sign')
%         toc
    end
    toc
    tic
    S2 = featuresign(B,X,lambda,S_init);
    toc
end