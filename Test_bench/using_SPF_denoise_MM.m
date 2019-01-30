function [Elmn_vec,d_recon] = using_SPF_denoise_MM(lambda_l,lambda_n,SNR,d_noisy,dlmn_GT,ms_XYZB)
%MM denoising algorithm for NCC distributed noise combined with radial and
%angular regularisation
%inputs:
% lambda_l  the angular LB regularisation parameter
% lambda_n  the radial regularisation parameter
% SNR       signal to noise ratio
% d_noisy   samples of the diffusion signal
% dlmn_GT   the GT SPF coefficients
% ms_XYZB   sampling scheme coordinates
%
%outputs:
% Elmn_vec          the vector of SPF coefficients - size (N+1)*(L+1)*(L+2)/2, aranged from n=0,...,N
% d_recon           the reconstructed diffusion signal

%author: Alice Bates
%date last altered: 13/02/18

 

%calculates the reconstruction accuracy for the proposed scheme for
%Gaussian mixture model
alpha=0.5;
N=3; % number of shells
L = [2 4 6 8]; %band-limit in each shell
Lmax = max(L);
Bmax = 4000e6;% in si units seconds per m
diff_time = 50e-3; 
qmax = sqrt(Bmax/(4*pi^2*diff_time));
x = roots(LaguerreGen(N+1,alpha)); % N+1 roots
%calculate scale factor
zeta = qmax^2/x(1);
x = flipud(x);

%roots in q-space
q_vec = sqrt(x*zeta);


%% calulate the coefficients
% MM Coefficient estimation
[Elmn_vec, d_recon,~, ~] = SPF_denoise_MM(d_noisy, SNR, lambda_l,lambda_n,L,N,x,zeta,dlmn_GT,ms_XYZB,diff_time);

  

% figure; plot(spectral_NRMSE_vec,'*');
% figure; plot(cost_func_vec,'o');


   

   



