function [ODF_4D] = LS_create_ODF_in_q_space(V,E)
%find the SPF basis matrix for samples V
% inputs:
% V         the sample locations to evaluate the ODF in q-space   
% E         the diffusion signal as a 4D matrix (x,y,z,image)
%
%outputs:
% ODF_4D    the ODF evaluated at V q-space coordinates for every voxel
%
%author: Alice Bates
%date last changed: 23/03/18

%regularization parameters, set based on SNR


lambda_l = 0.07;
lambda_n =0.04;
std = 1/27;%34.225;
%load('std_mat.mat');

%find b0 image
b0_img = mean(E(:,:,:,[1,134:142]),4);

alpha=0.5;
N=3; % number of shells
L = [2 4 6 8]; %band-limit in each shell
Lmax = max(L);
%Lmax=8; %for single shell ODF
num_SH = 0.5*(Lmax+1)*(Lmax+2);

alpha=0.5;
Bmax = 4000e6;% in si units seconds per m
diff_time = 50e-3; 
qmax = sqrt(Bmax/(4*pi^2*diff_time));
x = roots(LaguerreGen(N+1,alpha)); % N+1 roots
%calculate scale factor
zeta = qmax^2/x(1);
x = flipud(x);

%roots in q-space
q_vec = sqrt(x*zeta);

% %find real SH coeffs
V = [V(:,3),V(:,2),V(:,1)]; %[V(:,1),V(:,2),-V(:,3)];
%V_basis_mat = construct_SH_basis(Lmax, V, 2, 'real');
V_basis_mat = construct_SH_basis(Lmax, V, 2, 'complex');

%read in sampling coordinates
ms_XYZB =[];
for i=0:N
    load(['shell_',num2str(i+1),'.mat']);
    ms_XYZB = [ms_XYZB;XYZB];
end

    
S = construct_SPF_basis(Lmax, N, ms_XYZB, diff_time,zeta,0);
S = transpose(S);


%find the regularisation matrices
LB = [];
    for l=0:2:Lmax
        for m=-l:l
            LB = [LB; (l^2)*(l + 1)^2];
        end
    end
LB  = diag(repmat(LB,N+1,1));

n_vec = [];
for n = 0:N
    n_vec = [n_vec; repmat(n^2*(n+1)^2,num_SH,1)];
end
N_mat = diag(n_vec);



NUM_X_VOX = size(E,1);
NUM_Y_VOX = size(E,2);
NUM_Z_VOX = size(E,3);

%for every voxel, calculate the diffusion signal coefficients in the SPF
%basis Elmn. From these coefficients, find the ODF coefficients. Then find 
%the ODF at the points V by expanding the ODF SPF coefficients. 
ODF_4D = zeros(NUM_X_VOX,NUM_Y_VOX,NUM_Z_VOX,size(V,1));
for i=  45:70%65:82%%NUM_X_VOX
    for j= 72 %%NUM_Y_VOX
        for k = 53:70% 54:60 %% %NUM_Z_VOX
           % std = std_mat(i,j,k);
          
                y = E(i,j,k,2:133)/b0_img(i,j,k);
               y= reshape(y,132,1);
                
%                
% %                %% perform regularised LS
%                 SPFMatrixReg = (S'*S+lambda_l*LB+lambda_n*N_mat)\(S');
% 
%                 Elmn = SPFMatrixReg*y;
%                 
                
                 %% Or perform denoising
                % MM Rician Spherical Harmonic Coefficient estimation
                denoise_init = 1;
                [Elmn] = denoise_dwi_rician_mm(y, S, LB, N_mat, 1/std, lambda_l, lambda_n, denoise_init);
                
                
                
                %% put Elmn into matrix form
                Elmn_mat = zeros(N,num_SH);
                start = 1;
                for ind = 1:N+1
                   Elmn_mat(ind,1:0.5*(L(ind)+1)*(L(ind)+2)) = Elmn(start:start+0.5*(L(ind)+1)*(L(ind)+2)-1); 
                   start = start + start+0.5*(L(ind)+1)*(L(ind)+2);
                end
                    
               
                
                
                
                coeff_ODF = ODF_SH_coefficients(zeta, Elmn_mat,Lmax,N);
   
               
                
                
%                 
% %                 %% using ODF_CSA single shell
%                 load('shell_3.mat');
%                 y=E(i,j,k,23:67)/b0_img(i,j,k); %extract 3rd shell samples
%                 y = reshape(y,45,1);
%                 coeff_ODF = ODF_CSA_SS(lambda_l,y,XYZB);
                
                
                %% plot ODF
%                 ODF_XYZ = [V;-V];
%                 [basisV, ~,~] = construct_SH_basis (Lmax, ODF_XYZ, 2, 'real');
% 
%                 CSA_ODF_samples = basisV*coeff_ODF;

                %added min-max normalisation - not sure if need
                %min_CSA_ODF = min(CSA_ODF_samples); max_CSA_ODF = max(CSA_ODF_samples);
                %CSA_ODF_samples = (CSA_ODF_samples - min_CSA_ODF)/(max_CSA_ODF - min_CSA_ODF);

%                 %% - Plot ODF
%                 F=load('724_sphere_facets.txt');
%                 set(gcf, 'color', 'white');
%                 screen_size = get(0, 'ScreenSize');
%                 f1 = figure; hold on;
%                 set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%                 cameramenu;
%                 % -----
%                 Origen = [0 0 0.25]; % Voxel coordinates [x y z] where the ODF will be displayed
%                 angle = 0:.01:2*pi;
%                 R = ones(1,length(angle));
%                 polarm(angle,R,'k'); hold on;
%                 plot_ODF(CSA_ODF_samples,ODF_XYZ, F, Origen); 
%                 title('\fontsize{14} QBI_CSA', 'FontWeight','bold')
%                 axis equal; axis off;

                
                %find the q-space samples of the ODF
                 samples = real(V_basis_mat*coeff_ODF);
                 
                for ind = 1:length(samples)
                    if samples(ind) < 0
                        samples(ind) =0;
                    end
                end
                
                min_CSA_ODF = min(samples); %max_CSA_ODF = max(samples);
                samples = (samples - min_CSA_ODF);%/(max_CSA_ODF - min_CSA_ODF);
                
                ODF_4D(i,j,k,:) = samples;               
          
        end
    end
end

