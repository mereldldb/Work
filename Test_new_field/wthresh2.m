function res = wthresh(X,~,T)
    res = sign(X).*max(zeros(length(X),1),abs(X)-T);
end