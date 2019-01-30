function[Elmn_vec,d_recon_vec] = SPF_regularization(lambda_l, lambda_n,d_noisy)
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
%date last altered: 21/02/18



%calculates the reconstruction accuracy for the proposed scheme for
%Gaussian mixture model
alpha=0.5;
N=3; % radial degree, number of shells =N+1 
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
b_vec = q_vec.^2*4*pi^2*diff_time;


    
%% calulate the coefficients
[Elmn] = fSPFt(d_noisy,L,N,x,zeta,lambda_l,lambda_n);
d_recon = iSPFt(Elmn,L,N,x,zeta);
Elmn_vec = [];
for i=1:N+1
    Elmn_row = remove_odd_degree_SH(Elmn(i,:), Lmax);
    Elmn_vec = [Elmn_vec;transpose(Elmn_row)];
end

d_recon_vec = [];
 for i=1:N+1
    d_recon_vec= [d_recon_vec;transpose(d_recon(i,1:0.5*(L(i)+1)*(L(i)+2)))]; 
end




