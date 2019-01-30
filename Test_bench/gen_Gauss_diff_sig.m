function [y, y_noisy, noise_var] =  gen_Gauss_diff_sig(SNR, XYZB,varargin) %ANGLE,lambda_m,lambda_p
% gen_Gauss_diff_sig.m returns the diffusion signal sampled on XYZB with noise added.
% inputs:
% SNR       SNR on the B0 image 
% XYZB      The sampling scheme in Cartesian coordinates and the b-value as the 4th column. 
%optional inputs
% ANGLE     The crossing angle of the 2 fibres
% lambda_m  The diffusivity along the main tensor axis
% lambda_p  The diffusivity perpendicular to the main axis
%    
% outputs:
% y         the diffusion signal without noise
% y_noisy   the diffusion signal with noise added
% noise_var the variance of the noise = 1/SNR

    
    % create the voxel configuration
    VOXEL = MultiTensor();
    
% IF WANT 3 fibres    
%         VOXEL.M = 3;
%         VOXEL.f = [ 1/3, 1/3, 1/3 ];
%         VOXEL.lambda = [ 0.3 0.3 1.7 ; 0.3 0.3 1.7 ; 0.3 0.3 1.7 ]' * 1e-3;
%         
%          % create the rotation matrices to rotate the axis of the two fiber compartments
%         VOXEL.R(:,:,1) = VOXEL.ROTATION( 0, pi/2 ); %(PHI,THETA)
%         VOXEL.R(:,:,2) = VOXEL.ROTATION( pi/2, pi/2 ); %90 crossing angle
%         VOXEL.R(:,:,3) = VOXEL.ROTATION( 0, 0 );
    
    
    if nargin == 2
        % create two fiber compartments, with equal volume fraction and
        % diffusivities (default settings anyway)
        VOXEL.M = 2;
        VOXEL.f = [ 0.5, 0.5 ];
        VOXEL.lambda = [ 0.3 0.3 1.7 ; 0.3 0.3 1.7 ]' * 1e-3;
        
         % create the rotation matrices to rotate the axis of the two fiber compartments
        VOXEL.R(:,:,1) = VOXEL.ROTATION( 0, pi/2 ); %(PHI,THETA)
        VOXEL.R(:,:,2) = VOXEL.ROTATION( pi/2, pi/2 ); %90 crossing angle

 
    elseif nargin ==3
        ANGLE = varargin{1};
        % create two fiber compartments, with equal volume fraction and
        % diffusivities (default settings anyway)
        VOXEL.M = 2;
        VOXEL.f = [ 0.5, 0.5 ];
        VOXEL.lambda = [ 0.3 0.3 1.7 ; 0.3 0.3 1.7 ]' * 1e-3;
        
         % create the rotation matrices to rotate the axis of the two fiber compartments
        VOXEL.R(:,:,1) = VOXEL.ROTATION( 0, pi/2 ); %(PHI,THETA)
        VOXEL.R(:,:,2) = VOXEL.ROTATION( ANGLE, pi/2 );
            
    else %4 inputs
        lambda_m = varargin{1};
        lambda_p = varargin{2};
        % create two fiber compartments, with equal volume fraction and
        % diffusivities (default settings anyway)
        VOXEL.M = 1;
        VOXEL.f = 1;
        VOXEL.lambda = [ lambda_p lambda_p lambda_m ; lambda_p lambda_p lambda_m ]' * 1e-3;
    end

  
    % specify the signal-to-noise ratio
    sigma   = 1 / SNR; %as E0=1
    noise_var = sigma^2;

    [y, y_noisy ]= VOXEL.acquireWithScheme(XYZB,sigma);
    
