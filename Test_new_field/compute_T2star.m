function T2_inv = compute_T2star( Bz, gamma, npx )

%
% Construct T2_inv = gamma delta_Bz
% delta_Bz: change of magnetic field in pixel
%

D  = gallery('tridiag',npx,-0.5,0,0.5); % Difference quotient
I  = speye(npx,npx);
Nx = kron(I,D);
Ny = kron(D,I);

delta_Bz_x = Nx*Bz; % Jump in field strength in x-direction
delta_Bz_y = Ny*Bz; % Jump in field strength in y-direction
delta_Bz   = max(abs(delta_Bz_x),abs(delta_Bz_y)); % Take the maximum jump
T2_inv     = gamma*delta_Bz;

return
