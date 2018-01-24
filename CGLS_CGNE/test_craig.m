xc = A'*pcg(A*A',b,1e-15,1000);
%%
lambda = 1e-2;
Atil = [A lambda*eye(size(A,1))];
xs = pcg([A lambda*eye(size(A,1))]*[A lambda*eye(size(A,1))]',b,1e-15,1000);
xcraig = [A lambda*eye(size(A,1))]'*xs;
% imagesc(reshape(real(xcraig(1:4096)),64,64))
% colorbar

%% Wrong!!
lambda2 = 1e-4;
Atil2 = [A sqrt(lambda2)*eye(size(A,1)) zeros(size(A,1),size(F_hr,1)); ...
F_hr zeros(size(F_hr,1),size(A,1)) -eye(size(F_hr,1))];
btil = [b; zeros(size(F_hr,1),1)];
xv = pcg(Atil2*Atil2',btil,1e-20,1000);
xcraig2 = Atil2'*xv;
% figure(2)
% imagesc(reshape(real(xcraig2(1:4096)),64,64))
% colorbar
%%
subplot(2,2,1)
imagesc(reshape(real(xc),64,64))
title('Original Craig')
subplot(2,2,2)
imagesc(reshape(real(xcraig(1:4096)),64,64))
title('Saunders Craig')
subplot(2,2,3)
imagesc(reshape(real(xcraig2(1:4096)),64,64))
title('Extended Craig')
subplot(2,2,4)
imagesc(reshape(real(x_cgls2(:,end)),64,64))
title('CGLS')