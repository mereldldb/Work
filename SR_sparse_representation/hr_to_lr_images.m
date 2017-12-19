load('images')
for i = 1:10
    for m = 1:168
        for n = 1:168
            images_lr(n,m,i) = mean(mean(IMAGES((3*n-2:3*n),(3*m-2:3*m),i)));
        end
    end
end
save('images_lr','images_lr')