load('images.mat')
image_nr = 2;
X2 = IMAGES(:,:,image_nr);
patch_size_lr = 3;
patch_size_hr = 9;
L = patch_size_hr/patch_size_lr;
if L == 3
    load('images_lr.mat')
    x2 = images_lr(:,:,image_nr);
elseif L==2   
    load('images_lr_L2.mat')
    x2 = images_lr_L2(:,:,image_nr);
end

%X2 = X2(1:4*patch_size,1:4*patch_size);
npx = length(x2(1,:));
% lambda = .5;
% dict_size = 50;
% load('training_set')
% X_train = X(:,1:5000);
% D = train_dict(X_train,lambda,dict_size);
gamma = .01;
P = sparse(patch_size_lr^2,patch_size_lr^2);
matr_indices = reshape(1:npx^2,npx,npx);
counter = 1;

overlap = 0;
overlap_hr = 0;
beta = 0;
Xres = zeros(npx*L,npx*L);

for i = 1:patch_size_lr-overlap:npx-patch_size_lr+1
    for j = 1:patch_size_lr-overlap:npx-patch_size_lr+1
        indices_patch = matr_indices(j:j+patch_size_lr-1,i:i+patch_size_lr-1);
        if i == 1 && j == 1 %very first patch
            P = sparse(patch_size_hr^2,patch_size_hr^2);
        elseif j == 1 %first patch on the left
            P(overlap_hr+1:patch_size_hr,1:overlap_hr) = 0;
        elseif i == 1 %first patch from the top
            P(1:overlap_hr,overlap_hr+1:patch_size_hr) = 0;
        end
        indices_patch = reshape(indices_patch,1,[]);
        %F is for ordering patches as vectors
        F = sparse(patch_size_lr^2,npx^2);
        F(sub2ind([patch_size_lr^2,npx^2],1:patch_size_lr^2,indices_patch)) = 1;
        D_tilde = [D_lr; beta*P*D_hr];
        % for LR 3 and HR 9 we have 3*j-2:3*j+6,3*i-2:3*i+6
        y_tilde = [F*reshape(x2,[],1); beta*P*reshape(Xres(L*j-patch_size_lr+1:L*j+patch_size_hr-patch_size_lr,L*i-patch_size_lr+1:L*i+patch_size_hr-patch_size_lr),[],1)];
        alpha(:,counter) = admm(D_tilde,y_tilde,sparse(size(D_tilde,2),1),10,1e-10,20,1e-10,gamma,10*gamma);
%       alpha(:,counter) = ls_featuresign_sub(D_tilde,y_tilde,D_tilde'*D_tilde,D_tilde'*y_tilde,gamma);
        xres(j:j+patch_size_lr-1,i:i+patch_size_lr-1) = reshape(D_lr*alpha(:,counter),patch_size_lr,patch_size_lr);
        Xres(L*j-patch_size_lr+1:L*j+patch_size_hr-patch_size_lr,L*i-patch_size_lr+1:L*i+patch_size_hr-patch_size_lr) = reshape(D_hr*alpha(:,counter),patch_size_hr,patch_size_hr);        
        counter = counter + 1;
    end
end
%Hier nog backprojection toepassen, zie artikel
%%
figure(2)
subplot(2,2,1)
imagesc(X2)
colorbar
title('Model HR solution')
subplot(2,2,2)
imagesc(Xres)
colorbar
title('SR: HR solution')
subplot(2,2,3)
imagesc(x2)
colorbar
title('Model LR solution')
subplot(2,2,4)
imagesc(xres)
colorbar
title('SR: LR solution')