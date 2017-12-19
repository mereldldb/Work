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
sigma = 0.14e-9 * sqrt(bw_pulse); % standard deviation

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
tol_maj = 1e-6;
iter_maj = 10;
T = 0.05;

Dx = gallery('tridiag',npx,0,1,-1);
Ix = speye(npx);
Dy = gallery('tridiag',npy,0,1,-1);
Iy = speye(npy);
F = [kron(Iy,Dx);kron(Dy,Ix)]; 
F = F(sum(F,2)==0,:);
R = F'*F;
%R = speye(size(A,2));

lambda_t = 4e-17;%1e-15;%6e-16;%2e-16;

x_tv = admm_tv(A, b, sparse(npx*npy,1), speye(length(x_mod)),speye(length(b)), R, F, iter_admm, tol_admm, iter_CG, tol_CG, lambda_t, 10*lambda_t, x_mod);
%error2_tv = norm(max(0,real(x_tv))-x_mod,2);
%error1_tv = norm(max(0,real(x_tv))-x_mod,1);

% lambda_e = 8e-15;%3e-14;%9e-15;%4e-15; 
% order = 1;
% x_ep = majorization(A,b,speye(npx*npy),speye(length(b)),speye(length(b)), ...
%                            iter_maj,tol_maj,iter_fp,tol_fp,iter_CG,tol_CG,lambda_e,x_mod, sparse(npx*npy,1), npy, npx, T, order);
% error2_ep = norm(max(0,real(x_ep))-x_mod,2);
% error1_ep = norm(max(0,real(x_ep))-x_mod,1);

fig1 = figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
%subplot(1,2,1);
%imagesc(f, [0 1]); colorbar
%title('Model solution')

% subplot(2,2,1);
imagesc(reshape(real(x_tv),npy,npx),[0 1]); colorbar
axis square
% 
% subplot(2,2,2);
% imagesc(reshape(real(x_ep),npy,npx),[0 1]); colorbar
% axis square
