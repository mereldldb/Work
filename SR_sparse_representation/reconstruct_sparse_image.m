load('images.mat')
X2 = IMAGES(:,:,2);
patch_size = 8;
% X2 = X2(1:5*patch_size,1:5*patch_size);
npx = length(X2(1,:));
% lambda = .5;
% dict_size = 50;
% load('training_set')
% X_train = X(:,1:5000);
% D = train_dict(X_train,lambda,dict_size);
gamma = .2;
no_of_patches = (npx/patch_size)^2;
F = sparse(no_of_patches*patch_size^2,npx^2);
matr_indices = reshape(1:npx^2,npx,npx);
counter = 1;

for i = 1:patch_size:npx
    for j = 1:patch_size:npx
        indices = reshape(matr_indices(j:j+patch_size-1,i:i+patch_size-1),1,[]);
        %F is for ordering patches as vectors
        F(sub2ind([no_of_patches*patch_size^2,npx^2],counter:counter+patch_size^2-1,indices)) = 1;
        counter = counter + patch_size^2;
    end
end
alpha = featuresign(D,reshape(F*reshape(X2,[],1),patch_size*patch_size,no_of_patches),gamma);

Xtemp = D*alpha;
Xres = [];
% for p = 1:no_of_patches
%     
% end

% Xres = reshape(Xres, npx,npx)
counter2 = 0;
for k = 1:patch_size:npx
    for l = 1:patch_size:npx
        counter2 = counter2 + 1;
        Xres(l:l+patch_size-1,k:k+patch_size-1) = reshape(Xtemp(:,counter2),patch_size,patch_size);
    end
end
figure(1)
subplot(1,2,1)
imagesc(X2)
colorbar
subplot(1,2,2)
imagesc(Xres)
%imagesc([patch1 patch3; patch2 patch4])
colorbar