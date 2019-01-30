% Load dMRI data using SPM
path_to_data = 'process_real_data/opt_dim_multi.nii';
path_to_mask = 'process_real_data/opt_dim_multi_mask.nii';
path_to_save_data = 'process_real_data/ODF_folder';
Vdata = spm_vol(path_to_data);
E = spm_read_vols(Vdata);
MaskData = spm_vol(path_to_mask);
Mask = spm_read_vols(MaskData);

E = E.*Mask;
clear Mask MaskData;
% Load the spherical grid to evaluate the ODFs
V = load('sampling_and_reconstruction_schemes/On_the_sphere/724_shell.txt');
Vrec = V(1:362,:); % You can use this to save the data, what is Vrec? Half of the ODF samples? How different to ODF_4D?
  

%% add code to evaluate ODF in q-space. Need to find the ODF
%coefficients of the diffusion signal and evaluate them on the grid to make
%ODF_4D
ODF_4D = create_ODF_in_q_space(Vrec,E);

% Save ODFs using SPM
mkdir(path_to_save_data); % define path and create directory

Vaux = Vdata; % same as the diffusion data
Vaux = rmfield(Vaux,'private');
S = size(ODF_4D); % ODF_4D is the name of your ODFs in 4D
% ---
for i=1:S(4)
    VE(i) = Vaux(1); %just using to set up size 
    VE(i).fname = [path_to_save_data '/ODF_matrix_4D.nii'];
    VE(i).dt = [16 0];
    VE(i).n = [i 1];
    VE(i).dim = [S(1) S(2) S(3)];
    spm_write_vol(VE(i),squeeze(ODF_4D(:,:,:,i)));
end