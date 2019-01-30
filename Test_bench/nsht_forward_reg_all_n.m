function [Elm]=nsht_forward_reg_all_n(f,L,shell_num,lambda_l,lambda_n,n_vec)
% nsht_forward_reg_all_n.m - Computes the regularized forward spherical harmonic
% transform as a part of the SPF transform for every radial degree n. This
% is for a single-shell.
%
% Computes forward spherical harmonic transform based on the sampling
% scheme presented in paper:
%
%
% Inputs:
%
% f 		samples on the sphere
% L 		SH band-limit, max SH degree
% shell_num the number of the shell doing the transform over
% lambda_l  the angular LB regularisation parameter
% lambda_n	the radial regularisation parameter
% n_vec     the radial degrees to compute the regularised SHT over
%
%
% Outputs:
%
% Elm       a (N+1) x (L+1)^2 matrix containing the SH coefficients calculated for
% 			every radial order regularisation parameter.
%
%
%       
% Author: Alice Bates - adapted from code by Zubair Khalid
% Date last changed: 13/02/18
Elm = zeros(length(n_vec),(L+1)^2);

%read in Pm matrices for the current shell  
load(['Pm_shell_',num2str(shell_num),'.mat']); 

%perform SHT for each radial order n
for n_rad = n_vec
    [flm] = nsht_forward_reg_algo(f,L+1,lambda_l,lambda_n,n_rad,Pm_struct); %uses sig pro definition band-limit
    Elm(n_rad+1,:) = flm;
end


    