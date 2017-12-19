load('images.mat')
X2 = IMAGES(:,:,2);
patch_size = 8;
%X2 = X2(1:4*patch_size,1:4*patch_size);
npx = length(X2(1,:));
% lambda = .5;
% dict_size = 50;
% load('training_set')
% X_train = X(:,1:5000);
% D = train_dict(X_train,lambda,dict_size);
gamma = .2;
P = sparse(patch_size^2,patch_size^2);
matr_indices = reshape(1:npx^2,npx,npx);
counter = 1;

overlap = 2;
beta = 0;
Xres = zeros(npx,npx);

for i = 1:patch_size-overlap:npx-patch_size+1
    for j = 1:patch_size-overlap:npx-patch_size+1
        indices_patch = matr_indices(j:j+patch_size-1,i:i+patch_size-1);
        if i == 1 && j == 1 %very first patch
            P = sparse(patch_size^2,patch_size^2);
        elseif j == 1 %first patch on the left
            P(overlap+1:patch_size,1:overlap) = 0;
        elseif i == 1 %first patch from the top
            P(1:overlap,overlap+1:patch_size) = 0;
        end
        indices_patch = reshape(indices_patch,1,[]);
        %F is for ordering patches as vectors
        F = sparse(patch_size^2,npx^2);
        F(sub2ind([patch_size^2,npx^2],1:patch_size^2,indices_patch)) = 1;
        D_tilde = [D; beta*P*D];
        y_tilde = [F*reshape(X2,[],1); beta*P*reshape(Xres(indices_patch),[],1)];
        alpha(:,counter) = ls_featuresign_sub(D_tilde,y_tilde,D_tilde'*D_tilde,D_tilde'*y_tilde,gamma);
        Xres(j:j+patch_size-1,i:i+patch_size-1) = reshape(D*alpha(:,counter),patch_size,patch_size);        
        counter = counter + 1;
    end
end
figure(2)
subplot(1,2,1)
imagesc(X2)
colorbar
subplot(1,2,2)
imagesc(Xres)
colorbar