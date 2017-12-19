% function B = train_dict(X_train,lambda,dict_size)
load('sr_training_set')
hr_size = 9;
lr_size = 3;
dict_size = 500;
X_train = patches;
X_train = X_train(:,1:5000);
B = rand(length(X_train(:,1)),dict_size);
S = rand(dict_size,length(X_train(1,:)));
lambda = .6;
for i = 1:5
    disp('FS')
%     for g = 1:size(X_train,2)
%         S(:,g) = admm(B,X_train(:,g),sparse(size(B,2),1),10,1e-10,20,1e-10,lambda,10*lambda);
%     end
    S = featuresign(B,X_train,lambda,S);
    disp('Basis')
    B = learn_basis_from_author(X_train,S,B);
    err(i) = norm(X_train-B*S,'fro');
end
B_hr = hr_size*B(1:hr_size^2,:);
B_lr = lr_size*B(hr_size^2+1:end,:);