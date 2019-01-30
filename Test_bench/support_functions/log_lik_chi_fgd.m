function lik = log_lik_chi_fgd(y,c,B,Nc,noise_std2,N,Lmax,lambda_l, lambda_n)
% returns the value of the cost function for the iteration. The cost
% function is the Penalised maximum likelihood of the NCC distribution
% inputs:
% y 	 measurements of the diffusion signal
% c 	 coefficients of the diffusion signal
% B 	 basis matrix
% Nc  number of channels (Rician: Nc=1)
% noise_std2	 noise variance (Gaussian)
% N 		 Radial bandlimit for SPF basis
% max_L 	 SH bandlimit of the shell with largest SH bandlimit	
% lambda_l	 Laplace Beltrami angular regularization parameter
% lambda_n	 radial regularization parameter
%
% output
% lik 		the NLL cost function evaluated for the current iteration 
%
% author: Alice Bates - altered code from Divya Varadarajan

num_SH = 0.5*(Lmax+1)*(Lmax+2);

L_mat =[];
    for l=0:2:Lmax
        for m=-l:l
            L_mat = [L_mat; (l^2)*(l + 1)^2];
        end
    end
L_mat  = diag(repmat(L_mat,N+1,1));

n_vec = [];
for n = 0:N
    n_vec = [n_vec; repmat(n^2*(n+1)^2,num_SH,1)];
end
N_mat = diag(n_vec);

lik = real(sum((Nc-1)*log(B*c) +((B*c).^2)/(2*noise_std2)-(logbessi(Nc-1,(y.*(B*c))./noise_std2))) + lambda_n*transpose(c)*N_mat*c + lambda_l*transpose(c)*L_mat*c) ;
                       
