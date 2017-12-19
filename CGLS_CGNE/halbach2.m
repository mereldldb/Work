% Script to test MRI based on Halbach array
%
% Programmed by Melissa Wijchers, and Martin van Gijzen

%clear all; %close all;
close all

scrsz = get(0,'ScreenSize');

gamma = 267.513e6;     % rad/(sT)

%% Load magnetic field
load('field')

% Determine the frequency band for the field
BB_min  = min(min(BB)); freq_min = gamma*BB_min/(2*pi);
BB_max  = max(max(BB)); freq_max = gamma*BB_max/(2*pi);
fc_field = (freq_max+freq_min)/2;
bw_field = freq_max-freq_min;

disp(['Center frequency field = ',num2str(fc_field)]);
disp(['Bandwidth field = ',num2str(bw_field)]);

L = 8;
bw_pulse = 2000;
angles   = 0:10:350;   % rotation angles
FoV      = .12;%.04;%.02;
npx      = 64;
t_rec    = 1e-3; % Signal length
%% Menu
Title = 'MRI Parameters';
prompt = {'Pulse bandwidth (Hz):','Angles (degr):', ...
          'Field of View (m)', 'Number of pixels wide:','Recording time (s)'};
defaults = {num2str(bw_pulse),num2str(angles),num2str(FoV), ...
            num2str(npx), num2str(t_rec)};
lineNo=[1 40];
params=inputdlg(prompt,Title,lineNo,defaults);

bw_pulse = str2num(char(params(1)));
angles   = str2num(char(params(2)));
FoV      = str2num(char(params(3)));
npx      = str2num(char(params(4)));
t_rec    = str2num(char(params(5)));

npy = npx;

n_bands = ceil((bw_field)/bw_pulse);
f_min  = fc_field - n_bands*bw_pulse/2;
f_max  = fc_field + n_bands*bw_pulse/2;
% Compute center frequencies of the pulses:
fc_pulse = f_min+bw_pulse/2:bw_pulse:f_max-bw_pulse/2; 

disp(['Number of pulses = ',num2str(n_bands)]);

%% Generate phantom (by Merel)
f    =MRIphantom(npx);
x_mod=reshape(f,[],1);

%% Compute matrix
sigma = 1e-9;%0.14e-10 * sqrt(bw_pulse); % standard deviation

%sigma = 2e-9*sqrt(bw_pulse);
A = generate_matrix2(yy, zz, BB, gamma, fc_pulse, bw_pulse, angles, npx, FoV, t_rec );
[m,n] = size(A);
% Noise vector:
rng(0);
db = sigma*(randn(m,1)+sqrt(-1)*randn(m,1))/sqrt(2);

%%
disp('HR')
disp(['Number of equations  = ',num2str(m)]);
disp(['Number of pixels = ',num2str(n)]);

% Data vector:
d = A*x_mod;
b = d + db;

lambda = 1e-17;
tol_fp = 1e-100;
iter_fp = 20;
tol_CG = 1e-15;
iter_CG = 1000;
tol_admm = 1e-20;
iter_admm = 100;

Dx = gallery('tridiag',npx,0,1,-1);
Ix = speye(npx);
Dy = gallery('tridiag',npy,0,1,-1);
Iy = speye(npy);
F_hr = [kron(Iy,Dx);kron(Dy,Ix)];
% F = F(sum(F,2)==0,:);
%F_hr = kron(Iy,Dx)+kron(Dy,Ix);
R_hr = F_hr'*F_hr;

% R_hr = speye(length(x_mod));

lambda2 = 1e-15;
[x_cgls2, res1_cgls, res2_cgls, err2_cgls, obj2_cgls] = cgls(A,b,zeros(length(x_mod),1),speye(length(x_mod)),speye(length(b)),R_hr,iter_CG,tol_CG,lambda2,x_mod,1);
[x_cgne2, r_cgne2, res1_cgne, res2_cgne, err2_cgne, obj2_cgne] = cgne(A,b,zeros(length(b),1),speye(length(b)),speye(length(b)),R_hr,iter_CG,tol_CG,lambda2,x_mod,1);

%%
figure

% subplot(1,3,1);
% hold on;
% xlabel('Number of ITERATIONS')
% ylabel('log_{10} |r_k|')
% %title('Residual norm');
% grid on;
% plot(log10(res1_cgls),'b','LineWidth',2);
% plot(log10(res1_cgne),'r--','LineWidth',2);
% %plot(log10(rnrm_cgls2),'g','LineWidth',2);
% %plot(log10(rnrm_cgne2),'k--','LineWidth',2);
% legend('CGLS','CGNE');
% hold off;
% 
% subplot(1,3,2);
% hold on;
% xlabel('Number of ITERATIONS')
% ylabel('log_{10} |s_k|')
% %title('Residual norm normal 6equations');
% grid on;
% plot(log10(res2_cgls/res2_cgls(1)),'b','LineWidth',2);
% plot(log10(abs(res2_cgne/res2_cgne(1))),'r--','LineWidth',2);
% %plot(log10(snrm_cgls2/snrm_cgls2(1)),'g','LineWidth',2);
% %plot(log10(abs(snrm_cgne2{1}/snrm_cgne{1}(1))),'k--','LineWidth',2);
% legend('CGLS','CGNE');
% hold off;

subplot(1,2,1)
hold on 
xlabel('Number of ITERATIONS')
ylabel('Objective value')
%title('Relative Error');
grid on;
plot(obj2_cgls,'b','LineWidth',2);
plot(obj2_cgne,'r--','LineWidth',2);
%plot(log10(enrm_cgls2),'g','LineWidth',2);
%plot(log10(enrm_cgne2),'k--','LineWidth',2);
legend('CGLS','CGNE');
axis([0 50 0 1e-12])
hold off;

subplot(1,2,2);
hold on;
xlabel('Number of ITERATIONS')
ylabel('Error')
%title('Relative Error');
grid on;
plot(err2_cgls,'b','LineWidth',2);
plot(err2_cgne,'r--','LineWidth',2);
%plot(log10(enrm_cgls2),'g','LineWidth',2);
%plot(log10(enrm_cgne2),'k--','LineWidth',2);
legend('CGLS','CGNE');
axis([0 50 0 30])
hold off;

%%
[x_cgls, rnrm_cgls, snrm_cgls, enrm_cgls, obj_cgls] = cs_cgls(A,b,speye(length(x_mod)),speye(length(b)),F_hr,iter_fp,tol_fp,iter_CG,tol_CG,lambda,x_mod);
[x_cgne, rnrm_cgne, snrm_cgne, enrm_cgne, obj_cgne] = cs_cgne(A,b,speye(length(b)),speye(length(b)),F_hr,iter_fp,tol_fp,iter_CG,tol_CG,lambda,x_mod);

%b = d+db;

disp(['SNR = ',num2str(norm(d)/norm(db))]);

%%
figure

subplot(1,2,1)
hold on
xlabel('Number of ITERATIONS')
ylabel('log_{10}(Objective value)')
%title('Relative Error');
grid on;
%plot(log10(enrm_cgls2),'g','LineWidth',2);
%plot(log10(enrm_cgne2),'k--','LineWidth',2);
for k = 1:iter_fp
    plot((k-1)*iter_CG+1:k*iter_CG,log10(obj_cgls((k-1)*iter_CG+1:k*iter_CG)),'b','LineWidth',2);
    plot((k-1)*iter_CG+1:k*iter_CG,log10(obj_cgne((k-1)*iter_CG+1:k*iter_CG)),'r--','LineWidth',2);
end
y1=get(gca,'ylim');
for k = 1:iter_fp-1
    plot([k*iter_CG+.5 k*iter_CG+.5],y1,'k')
end
xlim([0 iter_CG*iter_fp])
legend('CGLS','CGNE');
hold off;

subplot(1,2,2);
hold on;
xlabel('Number of ITERATIONS')
ylabel('log_{10}(Error)')
%title('Relative Error');
grid on;
%plot(log10(enrm_cgls2),'g','LineWidth',2);
%plot(log10(enrm_cgne2),'k--','LineWidth',2);
for k = 1:iter_fp
    plot((k-1)*iter_CG+1:k*iter_CG,log10(enrm_cgls((k-1)*iter_CG+1:k*iter_CG)),'b','LineWidth',2);
    plot(((k-1)*iter_CG1:k*iter_CG,log10(enrm_cgne((k-1)*iter_CG+1:k*iter_CG)),'r--','LineWidth',2);
end
y1=get(gca,'ylim');
for k = 1:iter_fp-1
    plot([k*iter_CG+.5 k*iter_CG+.5],y1,'k')
end
xlim([0 iter_CG*iter_fp])
legend('CGLS','CGNE');
hold off;

%%
iter_CG2 = 100;

[x_admm, err_admm, obj_admm] = admm_tv(A,real(b),sparse(npx^2,1),speye(length(b)),speye(length(b)),R_hr,F_hr,iter_admm,tol_admm,iter_CG2,tol_CG,lambda,10*lambda,x_mod);
[x_admm2, err_admm2, obj_admm2] = admm_tv_2(A,real(b),sparse(length(b),1),speye(length(b)),speye(length(b)),R_hr,F_hr,iter_admm,tol_admm,iter_CG2,tol_CG,lambda,100*lambda,x_mod);

%%
figure

subplot(1,2,1)
hold on
xlabel('Number of ITERATIONS')
ylabel('log_{10}(Objective value)')
%title('Relative Error');
grid on;
%plot(log10(enrm_cgls2),'g','LineWidth',2);
%plot(log10(enrm_cgne2),'k--','LineWidth',2);
for k = 1:iter_admm
    plot((k-1)*iter_CG2+1:k*iter_CG2,log10(obj_admm((k-1)*iter_CG2+1:k*iter_CG2)),'b','LineWidth',2);
    plot((k-1)*iter_CG2+1:k*iter_CG2,log10(obj_admm2((k-1)*iter_CG2+1:k*iter_CG2)),'r--','LineWidth',2);
end
y1=get(gca,'ylim');
% for k = 1:iter_fp-1
%     plot([k*iter_CG+.5 k*iter_CG+.5],y1,'k')
% end
xlim([0 iter_CG2*iter_admm])
legend('CGLS','CGNE');
hold off;

subplot(1,2,2);
hold on;
xlabel('Number of ITERATIONS')
ylabel('log_{10}(Error)')
%title('Relative Error');
grid on;
%plot(log10(enrm_cgls2),'g','LineWidth',2);
%plot(log10(enrm_cgne2),'k--','LineWidth',2);
for k = 1:iter_admm
    plot((k-1)*iter_CG2+1:k*iter_CG2,log10(err_admm((k-1)*iter_CG2+1:k*iter_CG2)),'b','LineWidth',2);
    plot((k-1)*iter_CG2+1:k*iter_CG2,log10(err_admm2((k-1)*iter_CG2+1:k*iter_CG2)),'r--','LineWidth',2);
end
y1=get(gca,'ylim');
% for k = 1:iter_fp-1
%     plot([k*iter_CG+.5 k*iter_CG+.5],y1,'k')
% end
xlim([0 iter_CG2*iter_admm])
legend('CGLS','CGNE');
hold off;



% %%
% clear A b
% tol_CG = 1e-6;
% iter_CG = 100;
% tol_fp = 1e-6;
% iter_fp = 10;
% 
% Dx = gallery('tridiag',npx/L,0,1,-1);
% Ix = speye(npx/L);
% Dy = gallery('tridiag',npy/L,0,1,-1);
% Iy = speye(npy/L);
% F_lr = [kron(Iy,Dx);kron(Dy,Ix)]; 
% %F_lr = kron(Iy,Dx)+kron(Dy,Ix);
% %F = F(sum(F,2)==0,:);
% R_lr = F_lr'*F_lr;
% 
% lambda = 1e-15;
% %R_lr = speye(length(x_mod)/L^2);
% D = sparse(npx/L*npy/L,npx*npy);
%  
% for i = 1:npx/L
%     for j = 1:npy/L
%         X_y_indicesD = ((L*(i-1)+1):(L*(i-1)+L));
%         X_x_indicesD = (L*(j-1)+1):(L*(j-1)+L);
%  
%         [X_yD,X_xD] = meshgrid(X_y_indicesD,X_x_indicesD);
%         X_indicesD = X_xD+L*npy/L*(X_yD-1);
%         X_indicesD = reshape(X_indicesD,[],1);
%  
%         D_row = j+(i-1)*npy/L;
%         D(D_row,X_indicesD) = 1/L^2*ones(1,L^2);
%     end
% end
%  
% N = generate_matrix2(yy, zz, BB, gamma, fc_pulse, bw_pulse, angles, npx/L, FoV, t_rec );
% yy = [];
% for k = 1:9
%     y_mod{k} = D*reshape(imrotate(f,-angles(k),'bilinear','crop'),[],1);
%  
%     d = N*y_mod{k};
%  
% %     db_r = reshape(db,m/length(angles),length(angles));
% %     db_r = reshape(db_r(:,k:9:end-9+k),[],1);
%     b{k} = d+db;
%          
%     disp(['SNR = ',num2str(norm(d)/norm(db))]);
%  
%     load(sprintf('G_64pixels_%ddegrees',angles(k)))
%     %G{k} = generate_rotation_matrix(npx,npx,1,angles(k));
%     A((k-1)*npx^2/L^2+1:k*npx^2/L^2,:) = D*G;
%     
%     [y_cgls{k}, rnrm_y_cgls{k}, snrm_y_cgls{k}, enrm_y_cgls{k}] = cs_cgls(N,b{k},speye(length(y_mod{k})),speye(length(b{k})),F_lr,iter_fp,tol_fp,iter_CG,tol_CG,lambda,y_mod{k});
%     %y_cgls{k} = y_cgls{k}(:,end);
%     yy = [yy; y_cgls{k}];
% 
%     [y_cgne{k}, rnrm_y_cgne{k}, snrm_y_cgne{k}, enrm_y_cgne{k}] = cs_cgne(N,b{k},speye(length(b{k})),speye(length(b{k})),F_lr,iter_fp,tol_fp,iter_CG,tol_CG,lambda,y_mod{k});
% 
% end    
% 
% %[x_cgls, rnrm_cgls, snrm_cgls, enrm_cgls] = cgls(A,b,sparse(npx*npy,1),speye(length(x_mod)),speye(length(b)),R,iter_CG,tol_CG,lambda,x_mod,1);
% %[x_cgne, rnrm_cgne, snrm_cgne, enrm_cgne, ~] = cgne(A,b,zeros(length(b),1),speye(length(b)),speye(length(b)),R,iter_CG,tol_CG,lambda,x_mod,1);
% % %%
% % fig1 = figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
% % subplot(1,2,1)
% % imagesc(reshape(real(y_cgls{2}),npy/L,npx/L),[0 1]); colorbar
% % axis square
% % subplot(1,2,2)
% % imagesc(reshape(real(y_cgne{2}),npy/L,npx/L),[0 1]); colorbar
% % axis square
% 
% %%
% figure
% 
% subplot(1,3,1);
% hold on;
% xlabel('Number of ITERATIONS')
% ylabel('log_{10} |r_k|')
% %title('Residual norm');
% grid on;
% plot(log10(rnrm_y_cgls{2}),'b','LineWidth',2);
% plot(log10(rnrm_y_cgne{2}),'r--','LineWidth',2);
% %plot(log10(rnrm_cgls2),'g','LineWidth',2);
% %plot(log10(rnrm_cgne2),'k--','LineWidth',2);
% legend('CGLS','CGNE');
% hold off;
% 
% subplot(1,3,2);
% hold on;
% xlabel('Number of ITERATIONS')
% ylabel('log_{10} |s_k|')
% %title('Residual norm normal equations');
% grid on;
% plot(log10(snrm_y_cgls{2}/snrm_y_cgls{2}(1)),'b','LineWidth',2);
% plot(log10(abs(snrm_y_cgne{2}/snrm_y_cgne{2}(1))),'r--','LineWidth',2);
% %plot(log10(snrm_cgls2/snrm_cgls2(1)),'g','LineWidth',2);
% %plot(log10(abs(snrm_cgne2{1}/snrm_cgne{1}(1))),'k--','LineWidth',2);
% legend('CGLS','CGNE');
% hold off;
% 
% subplot(1,3,3);
% hold on;
% xlabel('Number of ITERATIONS')
% ylabel('log_{10} |x - x_k|')
% %title('Relative Error');
% grid on;
% plot(log10(enrm_y_cgls{2}),'b','LineWidth',2);
% plot(log10(enrm_y_cgne{2}),'r--','LineWidth',2);
% %plot(log10(enrm_cgls2),'g','LineWidth',2);
% %plot(log10(enrm_cgne2),'k--','LineWidth',2);
% legend('CGLS','CGNE');
% hold off;
% 
% %% SR2
% lambda = .00003;
% iter_fp = 10;
% iter_CG = 1000;
% tol_CG = 1e-10;
% tol_fp = 1e-10;
% tol_admm = 1e-6;
% iter_admm = 100;
% 
% [xx_cgls, rnrm_cgls, snrm_cgls, enrm_cgls] = cs_cgls(A,real(yy),speye(length(x_mod)),speye(length(yy)),F_hr,iter_fp,tol_fp,iter_CG,tol_CG,lambda,x_mod);
% [xx_cgne, rnrm_cgne, snrm_cgne, enrm_cgne] = cs_cgne(A,real(yy),speye(length(yy)),speye(length(yy)),F_hr,iter_fp,tol_fp,iter_CG,tol_CG,lambda,x_mod);
% xx_admm = admm_tv(A,real(yy),sparse(npx^2,1),speye(length(yy)),speye(length(yy)),R_hr,F_hr,iter_admm,tol_admm,iter_CG,tol_CG,lambda,10*lambda,x_mod);
% %%
% figure
% subplot(1,3,1)
% imagesc(reshape(xx_cgls,64,64))
% subplot(1,3,2)
% imagesc(reshape(xx_cgne,64,64))
% subplot(1,3,3)
% imagesc(reshape(xx_admm,64,64))
% 
% %%
% figure
% 
% subplot(1,3,1);
% hold on;
% xlabel('Number of ITERATIONS')
% ylabel('log_{10} |r_k|')
% %title('Residual norm');
% grid on;
% plot(log10(rnrm_cgls),'b','LineWidth',2);
% plot(log10(rnrm_cgne),'r--','LineWidth',2);
% %plot(log10(rnrm_cgls2),'g','LineWidth',2);
% %plot(log10(rnrm_cgne2),'k--','LineWidth',2);
% legend('CGLS','CGNE');
% hold off;
% 
% subplot(1,3,2);
% hold on;
% xlabel('Number of ITERATIONS')
% ylabel('log_{10} |s_k|')
% %title('Residual norm normal equations');
% grid on;
% plot(log10(snrm_cgls/snrm_cgls(1)),'b','LineWidth',2);
% plot(log10(abs(snrm_cgne/snrm_cgne(1))),'r--','LineWidth',2);
% %plot(log10(snrm_cgls2/snrm_cgls2(1)),'g','LineWidth',2);
% %plot(log10(abs(snrm_cgne2{1}/snrm_cgne{1}(1))),'k--','LineWidth',2);
% legend('CGLS','CGNE');
% hold off;
% 
% subplot(1,3,3);
% hold on;
% xlabel('Number of ITERATIONS')
% ylabel('log_{10} |x - x_k|')
% %title('Relative Error');
% grid on;
% plot(log10(enrm_cgls),'b','LineWidth',2);
% plot(log10(enrm_cgne),'r--','LineWidth',2);
% %plot(log10(enrm_cgls2),'g','LineWidth',2);
% %plot(log10(enrm_cgne2),'k--','LineWidth',2);
% legend('CGLS','CGNE');
% hold off;
