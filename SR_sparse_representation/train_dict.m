% function B = train_dict(X_train,lambda,dict_size)
load('training_set')
dict_size = 50;
X_train = X(:,1:5000);
%X_train = X(:,1:5000);
B = rand(length(X_train(:,1)),dict_size);
S = rand(dict_size,length(X_train(1,:)));
lambda = .5;
for i = 1:30
    S = featuresign(B,X_train,lambda,S);
    B = learn_basis(X_train,S,B);
    err(i) = norm(X_train-B*S,'fro');
end