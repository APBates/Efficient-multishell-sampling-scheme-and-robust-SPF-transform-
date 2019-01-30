function [dlmn, d_recon, spectral_NRMSE_vec, cost_func_vec] = SPF_denoise_MM(d_noisy, SNR, lambda_l,lambda_n,L,N,x,zeta,dlmn_GT, ms_XYZB,diff_time)
% SPF_denoise_MM.m iteratively computes the diffusion signal coefficients in the SPF basis while accounting for Rician or NCC distributed noise
%
% Inputs:
% d_noisy   Noisy DWI data for one voxel 
% SNR       =1/sigma, SNR of the DWI data
% lambda_l  angular LB regularization parameter
% lambda_n  radial regularization parameter
% L         vector of SH band-limits for each shell, max SH degree each shell
% N         radial band-limit SH basis, max radial degree
% x         location of shells in q^2/zeta
% zeta      scaling factor of SPF basis
% dlmn_GT   GT SPF coefficients
% ms_XYZB   sample locations 
% diff_time parameter of the sampling sequence
%
%% Outputs:
% dlmn              SPF coefficients vector
% d_recon           the reconstructed diffusion signal at sample locations
% spectral_NRMSE_vec the NRMSE in the SPF basis at each iteration
% cost_func_vec     the value of the ML function trying to be minimised at each iteration
%
% Author:  Alice Bates altered from code by Divya Varadarajan
% Date last altered: 13/02/18


Nc=1; % Number of channels
Lmax = max(L);
num_SH = 0.5*(Lmax+1)*(Lmax+2);

d_noisy_vec = [];
for i=1:(N+1)
    d_noisy_vec = [d_noisy_vec;transpose(d_noisy(i,1:0.5*(L(i)+1)*(L(i)+2)))];
end

noise_std = 1/SNR;
noise_std2 = noise_std^2;
dovernoise2 = d_noisy_vec/noise_std2;

%generate imaginary SPF basis
S = construct_SPF_basis(Lmax, N, ms_XYZB, diff_time,zeta,0);
S = transpose(S);
 
% initalisation
dnew_init = fSPFt(d_noisy,L,N,x,zeta,lambda_l,lambda_n);
%remove odd degree SH and turn into a vector.
Elmn_vec = [];
for i=1:(N+1)
    Elmn_row = remove_odd_degree_SH(dnew_init(i,:), Lmax);
    Elmn_vec = [Elmn_vec;transpose(Elmn_row)];
end

dnew_init = Elmn_vec;
clear Elmn_vec;


% Coeff initialization
dnew = dnew_init;

% MM iteration
grad_at_xo = @(xo,dovernoise2) noise_std2*dovernoise2.*real(exp(logbessi(Nc,xo.*dovernoise2)...
                   -logbessi(Nc-1,xo.*dovernoise2)));


spectral_NRMSE_vec = [];
NRMSE = norm(dnew-dlmn_GT)/norm(dlmn_GT);
spectral_NRMSE_vec = [spectral_NRMSE_vec, NRMSE];
cost_func_vec =[];


neg_log_lik = log_lik_chi_fgd(d_noisy_vec, dnew,S,Nc,noise_std2,N,Lmax,lambda_l, lambda_n);
cost_func_vec = [cost_func_vec; neg_log_lik];

do = dnew;
iter3 = 0;
while(1) %until convergence
    iter3 = iter3+1
    g= grad_at_xo(S*do,dovernoise2);
    %put updated samples into matrix as fSPFt needs them in this form
    g_mat = zeros((N+1),num_SH);
    start=1;
    for i=1:(N+1)
       g_mat(i,1:0.5*(L(i)+1)*(L(i)+2)) = g(start:start+0.5*(L(i)+1)*(L(i)+2)-1);
       start = start + 0.5*(L(i)+1)*(L(i)+2);
    end
    
    dnew= fSPFt(g_mat,L,N,x,zeta,lambda_l,lambda_n);
    %remove odd degree SH and turn into a vector
    Elmn_vec = [];
    for i=1:(N+1)
        Elmn_row = remove_odd_degree_SH(dnew(i,:), Lmax);
        Elmn_vec = [Elmn_vec;transpose(Elmn_row)];
    end

    dnew = Elmn_vec;
    clear Elmn_vec;
    
    change = norm(do-dnew)/norm(do);
 
    
    %work out the cost function - should be decreasing
    clear lik;
    neg_log_lik = log_lik_chi_fgd(d_noisy_vec,dnew,S,Nc,noise_std2,N,Lmax,lambda_l, lambda_n);
    neg_log_lik = sum(neg_log_lik);
    cost_func_vec = [cost_func_vec; neg_log_lik];
    
    do = dnew;
    NRMSE = norm(dnew-dlmn_GT)/norm(dlmn_GT);
    spectral_NRMSE_vec = [spectral_NRMSE_vec; NRMSE];
    
    if (change < 1e-3 || iter3 >50) %when converged
        disp(['Coming out of while loop, iterations -' num2str(iter3)]);
        numel(find(sum(S*dnew<0,1)))
        break;
    end
    
   
end;

dlmn = dnew;
d_recon = S*dnew;
end

