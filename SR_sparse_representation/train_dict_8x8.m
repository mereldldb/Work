% function B = train_dict(X_train,lambda,dict_size)
load('training_set_3x3')
dict_size = 50;
X_train = patches;
%X_train = X(:,1:5000);
B = rand(length(X_train(:,1)),dict_size);
S = rand(dict_size,length(X_train(1,:)));
lambda = 1;
for i = 1:1
    S = featuresign(B,X_train,lambda,S);
    B = learn_basis(X_train,S,B);
    err(i) = norm(X_train-B*S,'fro');
end