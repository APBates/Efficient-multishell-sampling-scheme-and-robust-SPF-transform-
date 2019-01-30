function[Elmn_mat,d_recon_vec] = SPF_regularization_real_data(lambda_l, lambda_n,d_noisy_vec)
%Calulates the regularized recontruction of the diffusion signal d_recon and the SPF coefficients Elmn_vec using the
%spectral multi-shell sampling scheme
%inputs:
% lambda_l	the angular LB regularisation parameter
% lambda_n	the radial regularisation parameter
% d_noisy   Diffusion signal measurements
%
%outputs:
% Elmn_vec			the vector of SPF coefficients - size (N+1)*(L+1)*(L+2)/2, aranged from n=0,...,N
% d_recon_vec       the reconstructed diffusion signal 
% locations
%author: Alice Bates
%date last altered: 06/03/18



%calculates the reconstruction accuracy for the proposed scheme for
%Gaussian mixture model
alpha=0.5;
N=3; % radial degree, number of shells =N+1 
L = [2 4 6 8]; %band-limit in each shell
Lmax = max(L);
num_SH = 0.5*(Lmax+1)*(Lmax+2);
Bmax = 4000e6;% in si units seconds per m
diff_time = 50e-3; 
qmax = sqrt(Bmax/(4*pi^2*diff_time));
x = roots(LaguerreGen(N+1,alpha)); % N+1 roots
%calculate scale factor
zeta = qmax^2/x(1);
x = flipud(x);

%roots in q-space
q_vec = sqrt(x*zeta);
b_vec = q_vec.^2*4*pi^2*diff_time;

d_noisy_mat = zeros((N+1),num_SH);
start=1;
for i=1:(N+1)
    d_noisy_mat(i,1:0.5*(L(i)+1)*(L(i)+2)) = d_noisy_vec(start:start+0.5*(L(i)+1)*(L(i)+2)-1);
    start = start + 0.5*(L(i)+1)*(L(i)+2);
end
    
%% calulate the coefficients
[Elmn] = fSPFt(d_noisy_mat,L,N,x,zeta,lambda_l,lambda_n);
%d_recon = iSPFt(Elmn,L,N,x,zeta);

 d_recon_vec = [];
%  for i=1:N+1
%     d_recon_vec= [d_recon_vec;transpose(d_recon(i,1:0.5*(L(i)+1)*(L(i)+2)))]; 
%  end
 
%convert to real SH
Elmn_mat = zeros(N, num_SH); 
for i=1:N+1
    Elmn_row = remove_odd_degree_SH(Elmn(i,1:(L(i)+1)^2),L(i));
    Elmn_mat(i,1:0.5*(L(i)+1)*(L(i)+2)) = Elmn_row; %complex_to_real_SH(Elmn_row, L(i));
end



