function [ODF_4D] = SS_create_ODF_in_q_space(V,E)
%find the SPF basis matrix for samples V
% inputs:
% V         the sample locations to evaluate the ODF in q-space   
% E         the diffusion signal as a 4D matrix (x,y,z,image)
%
%outputs:
% ODF_4D    the ODF evaluated at V q-space coordinates for every voxel
%
%author: Alice Bates
%date last changed: 23/02/18

%regularization parameters, set based on SNR


lambda_l = 0.004;
std = 1/27;%34.225;
%load('std_mat.mat');
%load('optimal_reg.mat');

%find b0 image
b0_img = mean(E(:,:,:,[1,134:142]),4);


Lmax = 8;


b = 4000e6;% in si units seconds per m

%calculate scale factor
zeta = qmax^2/x(1);
x = flipud(x);

% %find real SH coeffs
V = [V(:,3),V(:,2),V(:,1)]; %[V(:,1),V(:,2),-V(:,3)];
%V_basis_mat = construct_SH_basis(Lmax, V, 2, 'real');
V_basis_mat = construct_SH_basis(Lmax, V, 2, 'complex');


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

                 load('shell_3.mat');
                 y=E(i,j,k,23:67)/b0_img(i,j,k); %extract 3rd shell samples
                 y = reshape(y,45,1);
                 coeff_ODF = ODF_CSA_SS_opt_dim_coeffs(y,lambda_l,XYZB,std);

                
                %find the q-space samples of the ODF
                 samples = real(V_basis_mat*coeff_ODF);
                 
                for ind = 1:length(samples)
                    if samples(ind) < 0
                        samples(ind) =0;
                    end
                end

                
                ODF_4D(i,j,k,:) = samples;               
          
        end
    end
end

