function [A,T2inv] = generate_matrix2( x, y, Bz, gamma, freq, bw, angles, npx, FoV, t_rec )

%
% Input:
%    Bz:     magnetic field
%    gamma:  gyromagnetic ratio [rad/Ts]
%    freq:   center frequencies of the pulses (can be multiple)
%    bw:     bandwidth of the pulses (the same for all)
%    angles: rotation angles of the magnetic field
%    npx:    number of pixels in one direction (phantom is square)
%    FoV:    Field of view (in one direction)
%    t_rec:  recording time
%

fs    = bw;                % Sample frequency for basebanded signal 
N     = round(t_rec*fs);   % number of timesteps
dt    = t_rec/N;
t     = [0:dt:(N-1)*dt]';

%% interpolated field Bz (value per pixel)
sizefield=size(Bz);
if length(FoV) == 1
    FoV = [FoV FoV];
end
dx=FoV(1)/(sizefield(1)-1); dy=FoV(2)/(sizefield(2)-1); dz = 0.005; % 5mm
%X = [0:dx:FoV(1)]'; Y = [0:dy:FoV(2)]';
X = x';
Y = y;
[X,Y] = meshgrid(X,Y);
dx_int = FoV(1)/(npx-1); dy_int = FoV(2)/(npx-1);
Xq = [-FoV(1)/2:dx_int:FoV(1)/2]'; Yq = [-FoV(2)/2:dy_int:FoV(2)/2];
[Xq,Yq] = meshgrid(Xq,Yq);
Bz_int = interp2(X,Y,Bz,Xq,Yq,'cubic');

% imagesc(100*Xq(1,:),100*Yq(:,1),1000*Bz_int)
% colorbar
% xlabel('x (cm)')
% ylabel('y (cm)')
% title('B (mT)')
% axis square

n_angles = length(angles);
n_freq   = length(freq);

n = npx*npx;

indx = zeros(N*n_angles*n,1);
indy = zeros(N*n_angles*n,1);
A    = zeros(N*n_angles*n,1);

m  = 0;  % Number of equations
nz = 0;  % Number of nonzeros

tic
for i=1:n_angles

% Generate the rotated magnetic field
   angle  = angles(i);
   Bz_rot = imrotate(Bz_int,angle,'bilinear','crop');
   Bz_rot = reshape(Bz_rot,n,1);
   omega0 = gamma *Bz_rot;
   M0     = 0.0031*Bz_rot; % Magnetization

%% Construct T2 
   T2inv = compute_T2star( Bz_rot, gamma, npx );

   for i_freq = 1:n_freq

% Generate the new equations

% Get center frequency
      fc = freq(i_freq);

% c: Receiver sensitivity times pixel volume
% Based on rough estimate in note by Andrew:
      c = (60e-6)*(dx_int*dy_int*dz);

%% Determine which pixels are excited by pulse
      ind = find( (omega0 < 2*pi*(fc+bw/2)) & (omega0 > 2*pi*(fc-bw/2)));

% B0: magnetic field strength that corresponds to the center frequency
      B0 = 2*pi*fc/gamma;

% dBz: basebanded magnetic field 
      dBz=Bz_rot-B0;

%% Compute new block Ab of A
%%    A_new = c*omega*exp(-t/T2*)*exp(-i*gamma*dBz*t)

% Factor c*omega*M0 (diagonal matrix, scaling of the unknowns):
      CWM0 = c*sparse( ind, ind, omega0(ind).*M0(ind), n, n );

% Factor exp(-t*i*gamma*dBz): Spins
% Factor exp(-t/T2): Damping of the signal
% Together: exp(-t(i*gamma*dBz+1/T2)):
% For efficiency first compute (1i*gamma)*dBz+T2inv and make it sparse

      z        = (1i*gamma)*dBz+T2inv;
      Z        = sparse( 1, ind, z(ind), 1, n );
      A_new    = exp(-t*Z)*CWM0;

% Put the nonzero coefficients and corresponding indices in A
      [ix,iy,a] = find(A_new);
      indx(nz+1:nz+length(ix)) = ix+m;
      indy(nz+1:nz+length(iy)) = iy;
      A(nz+1:nz+length(a))     = a;

      m = m + N;
      nz  = nz + length(ix);
   end
end
toc

% Make A sparse matrix
A = sparse(indx(1:nz),indy(1:nz),A(1:nz),m,n);

% Throw away all the zero rows:
sumA = sum(abs(A),2);
maxA = max(sum(A));
ind = sumA/maxA > 1e-8;
%A = A(ind,:);
return