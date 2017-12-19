% Script to test MRI based on Halbach array
%
% Programmed by Melissa Wijchers, and Martin van Gijzen

%clear all; %close all;
close all

scrsz = get(0,'ScreenSize');

gamma = 267.513e6;     % rad/(sT)

%% Load magnetic field
%load('Bz');
load('field')

% Determine the frequency band for the field
BB_min  = min(min(BB)); freq_min = gamma*BB_min/(2*pi);
BB_max  = max(max(BB)); freq_max = gamma*BB_max/(2*pi);
fc_field = (freq_max+freq_min)/2;
bw_field = freq_max-freq_min;

disp(['Center frequency field = ',num2str(fc_field)]);
disp(['Bandwidth field = ',num2str(bw_field)]);

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
L = 1;
%sigma = 0.14e-9 * sqrt(bw_pulse); % standard deviation
sigma = 1.4e-9 * sqrt(bw_pulse); % standard deviation

A = generate_matrix2(yy, zz, BB, gamma, fc_pulse, bw_pulse, angles, npx, FoV, t_rec );
[m,n] = size(A);
% Noise vector:
db = sigma*(randn(m,1)+sqrt(-1)*randn(m,1))/sqrt(2);

%%
disp('HR')
disp(['Number of equations  = ',num2str(m)]);
disp(['Number of pixels = ',num2str(n)]);

% Data vector:
d = A*x_mod;

b = d+db;

disp(['SNR = ',num2str(norm(d)/norm(db))]);
%%
tol_CG = 1e-6;
iter_CG = 10;
tol_admm = 1e-6;
iter_admm = 10;
tol_fp = 1e-6;
iter_fp = 10;

Dx = gallery('tridiag',npx,0,1,-1);
Ix = speye(npx);
Dy = gallery('tridiag',npy,0,1,-1);
Iy = speye(npy);
F_hr = [kron(Iy,Dx);kron(Dy,Ix)]; 
F_hr = F_hr(sum(F_hr,2)==0,:);
R_hr = F_hr'*F_hr;

lambda_hr = 1e-15;%1e-15;%6e-16;%2e-16;

x_tv = admm_tv(A, b, sparse(npx*npy,1), speye(length(x_mod)),speye(length(b)), R_hr, F_hr, iter_admm, tol_admm, iter_CG, tol_CG, lambda_hr, 10*lambda_t, x_mod);
error2_tv = norm(max(0,real(x_tv))-x_mod,2);
error1_tv = norm(max(0,real(x_tv))-x_mod,1);

% lambda_e = 8e-15;%3e-14;%9e-15;%4e-15; 
% order = 1;
% x_ep = majorization(A,b,speye(npx*npy),speye(length(b)),speye(length(b)), ...
%                            iter_maj,tol_maj,iter_fp,tol_fp,iter_CG,tol_CG,lambda_e,x_mod, sparse(npx*npy,1), npy, npx, T, order);
% error2_ep = norm(max(0,real(x_ep))-x_mod,2);
% error1_ep = norm(max(0,real(x_ep))-x_mod,1);

fig1 = figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
subplot(2,2,1);
imagesc(f, [0 1]); colorbar
title('Model solution')
axis square

 subplot(2,2,2);
imagesc(reshape(real(x_tv),npy,npx),[0 1]); colorbar
axis square
% 
% subplot(2,2,2);
% imagesc(reshape(real(x_ep),npy,npx),[0 1]); colorbar
% axis square
%%
disp('SR 4 measurements')
clear A b
L = 8;
D = sparse(npx/L*npy/L,npx*npy);
 
for i = 1:npx/L
    for j = 1:npy/L
        X_y_indicesD = ((L*(i-1)+1):(L*(i-1)+L));
        X_x_indicesD = (L*(j-1)+1):(L*(j-1)+L);
 
        [X_yD,X_xD] = meshgrid(X_y_indicesD,X_x_indicesD);
        X_indicesD = X_xD+L*npy/L*(X_yD-1);
        X_indicesD = reshape(X_indicesD,[],1);
 
        D_row = j+(i-1)*npy/L;
        D(D_row,X_indicesD) = 1/L^2*ones(1,L^2);
    end
end

N = generate_matrix2(yy, zz, BB, gamma, fc_pulse, bw_pulse, angles(1:9:28), npx/L, FoV, t_rec );
 
for k = 1:9
    y_mod{k} = D*reshape(imrotate(f,-angles(k),'bilinear','crop'),[],1);
 
    d = N*y_mod{k};
 
    db_r = reshape(db,m/length(angles),length(angles));
    db_r = reshape(db_r(:,k:9:end-9+k),[],1);
    b{k} = d+db_r;
         
    disp(['SNR = ',num2str(norm(d)/norm(db_r))]);
 
    load(sprintf('G_64pixels_%ddegrees',angles(k)));
    A((k-1)*npx^2/L^2+1:k*npx^2/L^2,:) = D*G;
end

lambda_4 = 0;

Dx = gallery('tridiag',npx/L,0,1,-1);
Ix = speye(npx/L);
Dy = gallery('tridiag',npy/L,0,1,-1);
Iy = speye(npy/L);
F_lr = [kron(Iy,Dx);kron(Dy,Ix)]; 
F_lr = F_lr(sum(F_lr,2)==0,:);
R_lr = F_lr'*F_lr;

yt = [];

for l = 1:9
      %[y{l},err] = admm_tv(N, b{l}, sparse(npx*npy/L^2,1), speye(npx^2/L^2),speye(length(b{l})), R_lr, F_lr, iter_admm, tol_admm, iter_CG, tol_CG, lambda, 10*lambda, y_mod{l});
 
      y{l} = cgls(N, b{l}, sparse(npx*npy/L^2,1), speye(npx^2/L^2),speye(length(b{l})), R_lr, iter_CG, tol_CG, lambda_4, y_mod{l},1);
     
      y{l} = y{l}(:,end);
      error2(l) = norm(real(y{l})-y_mod{l},2);
      error1(l) = norm(real(y{l})-y_mod{l},1);
 
%       subplot(3,3,l)
%       axis square
%       imagesc(reshape(real(y{l}),npx/L,npx/L),[0 1]); colorbar
 
      yt = [yt; y{l}];
end
 
bb = real(yt);
bb = max(0,bb);

%%
lambda_sr4 = .003;

x4_sr_tv = admm_tv(A, bb, sparse(npx*npy,1), speye(length(x_mod)),speye(length(bb)), R_hr, F_hr, iter_admm, tol_admm, iter_CG, tol_CG, lambda_sr4, 10*lambda_sr4, x_mod);

error2_sr_tv = norm(max(0,real(x4_sr_tv))-x_mod,2);
error1_sr_tv = norm(max(0,real(x4_sr_tv))-x_mod,1);
%     [error2_sr_tv, error1_sr_tv]

figure(1)
subplot(2,2,3)
imagesc(reshape(real(x4_sr_tv),npy,npx),[0 1]); colorbar
axis square


%%
disp('SR 9 measurements')

N = generate_matrix2(yy, zz, BB, gamma, fc_pulse, bw_pulse, angles, npx/L, FoV, t_rec );
 
for k = 1:9
    y_mod{k} = D*reshape(imrotate(f,-angles(k),'bilinear','crop'),[],1);
 
    d = N*y_mod{k};
    b{k} = d + db;
    %db_r = reshape(db,m/length(angles),length(angles));
    %db_r = reshape(db_r(:,k:9:end-9+k),[],1);
    %b{k} = d+db_r;
         
    disp(['SNR = ',num2str(norm(d)/norm(db))]);
 
    %G = generate_rotation_matrix(npx,npx,1,angles(k));
    load(sprintf('G_64pixels_%ddegrees',angles(k)));
    A((k-1)*npx^2/L^2+1:k*npx^2/L^2,:) = D*G;
end
%%
fig1 = figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
     
lambda_9 = 0;
yt = [];
    
for l = 1:9
      %[y{l},err] = admm_tv(N, b{l}, sparse(npx*npy/L^2,1), speye(npx^2/L^2),speye(length(b{l})), R_lr, F_lr, iter_admm, tol_admm, iter_CG, tol_CG, lambda_tv, 10*lambda_tv, y_mod{l});
      
      y{l} = cgls(N, b{l}, sparse(npx*npy/L^2,1), speye(npx^2/L^2),speye(length(b{l})), R_lr, iter_CG, tol_CG, lambda_9, y_mod{l}, 1);
      y{l} = y{l}(:,end);
      
      error2(l) = norm(real(y{l})-y_mod{l},2);
      error1(l) = norm(real(y{l})-y_mod{l},1);
 
      subplot(3,3,l)
      axis square
      imagesc(reshape(real(y{l}),npx/L,npx/L),[0 1]); colorbar
 
      yt = [yt; y{l}];
end
 
bb = real(yt);
bb = max(0,bb);
 
%%
Dx = gallery('tridiag',npx,0,1,-1);
Ix = speye(npx);
Dy = gallery('tridiag',npy,0,1,-1);
Iy = speye(npy);
F = [kron(Iy,Dx);kron(Dy,Ix)]; 
F = F(sum(F,2)==0,:);
R = F'*F;
 
lambda_sr9 = .001;
% contin = 1;
 
% while ( contin ) 
%     Title = 'SR solver parameters';
%     prompt = {'Total variation lambda:', 'Edge-preserving lambda:'};
%     defaults = {num2str(lambda_tv), num2str(lambda_ep)};
%     lineNo=[1 40];
%     params=inputdlg(prompt,Title,lineNo,defaults);
%     if ( isempty(params) ) 
%        contin = 0; 
%        break; 
%     end
%  
%     lambda_tv = str2num(char(params(1)));
%     lambda_ep = str2num(char(params(2)));
 
    x9_sr_tv = admm_tv(A, bb, sparse(npx*npy,1), speye(length(x_mod)),speye(length(bb)), R_hr, F_hr, iter_admm, tol_admm, iter_CG, tol_CG, lambda_sr9, 10*lambda_sr9, x_mod);
 
    error2_sr_tv = norm(max(0,real(x9_sr_tv))-x_mod,2);
    error1_sr_tv = norm(max(0,real(x9_sr_tv))-x_mod,1);
%     [error2_sr_tv, error1_sr_tv]
 
    figure(1)
    subplot(2,2,4)
    imagesc(reshape(real(x9_sr_tv),npy,npx),[0 1]); colorbar
    axis square
 

% end
