load('images')
lr_size = 2;
hr_size = 4;
patches = zeros(lr_size^2+hr_size^2,50000);
for p = 1:50000
    k = unidrnd(9) + 1;
    i = unidrnd(511-hr_size) + 1;
    j = unidrnd(511-hr_size) + 1;
    patch_hr = IMAGES(i:i+hr_size-1,j:j+hr_size-1,k);
    counter = 1;
    for m = 1:lr_size:hr_size
        for n = 1:lr_size:hr_size
            patch_lr(counter,1) = mean(mean(patch_hr(n:n+lr_size-1,m:m+lr_size-1)));
            counter = counter+1;
        end
    end
    patch_hr = reshape(patch_hr,[],1);
    patches(:,p) = [1/hr_size*patch_hr; 1/lr_size*patch_lr];
end
save('sr_training_set_2and4','patches')