function [flm] = nsht_forward_reg_algo(f,L_sig,lambda_l,lambda_n,n_rad,Pm_struct)
% nsht_forward_reg_algo - Computes the regularized forward spherical harmonic
% transform as a part of the SPF transform.
%
% Computes forward spherical harmonic transform based on the sampling
% scheme presented in paper:
%
%
% Inputs:
%
% Outputs:
%
% where L_sig is the harmonic band-limit used by the signal processing community, f is the vector of L_sig^2 evaluated over the sampling scheme - this is the sig pro definition of band-limit = Lmax+1 
% and flm is the vector of L^2 harmonic coefficients
%
%
%       
% Author: Alice Bates - adapted from code by Zubair Khalid

if ~isreal(L_sig)
      error('Harmonic band-limit must be real');
end

if ~(sum(size(f))== length(f)+1) 
      error('f must be a vector and not a matrix');
end


    flm = zeros(1,(L_sig)^2);

    [~, FI] = nsht_sampling_points(L_sig); 

    fft_v = zeros(1,(L_sig)^2);

    for m=L_sig-1:-2:0
        P_mat_1 = Pm_struct{m+1}; 

        if m~=0
            P_mat_2=Pm_struct{m};
         end

        %uniform sampling - FFT version
        fft_v_bit = fft(f((m^2-m+2)/2:(m+1)*(m+2)/2))/length(f((m^2-m+2)/2:(m+1)*(m+2)/2));      

        fft_v(m^2+1:(m+1)^2) = fft_v_bit;

        flm_t_1 = nsht_inversion( P_mat_1, fft_v,m,L_sig,n_rad,lambda_l,lambda_n);
        if m~=0
            flm_t_2 = nsht_inversion( P_mat_2, fft_v,m-1,L_sig,n_rad,lambda_l,lambda_n); %coefficients for m-1
        end

        % perform spatial elimitation
        [ f ] = nsht_spatial_elimination(P_mat_1, f(1:m*(m-1)/2),flm_t_1, FI(1:m*(m-1)/2), m );
         if m~=0
            [ f ] = nsht_spatial_elimination(P_mat_2, f(1:m*(m-1)/2),flm_t_2, FI(1:m*(m-1)/2), m-1 );
         end

        %store coefficients calculated
        for el=m:1:L_sig-1
           flm(el^2+el+m+1) = flm_t_1(el-m+1);
           flm(el^2+el-m+1) = (-1)^m*conj(flm_t_1(el-m+1));
        end

       if m~=0 
           n = m-1;
           for el=n:1:L_sig-1

               flm(el^2+el+n+1) = flm_t_2(el-n+1);
               flm(el^2+el-n+1) = (-1)^n*conj(flm_t_2(el-n+1));
           end
       end

    end
    
end


function f_t = nsht_spatial_elimination(P_mat, f_t, flm_t,  FI_t, m)
% nsht_spatial_elimination - removes the m and -ve m order coefficients from the signal
% when flm is passed, assuming fl(-m) are for real signal. Input to
% function would be spatial signal 2m+1 coefficents and order m, TT and FF. The output
% would be spatial signal with coefficients removed.

% f_t = truncated signal 
% flm_t spectrum of order m contains L-m+1 coefficients, the first one
% corresponds to degree el=m and the last one for degree el=L-1

% THETA_t, FI_t, the truncated spatial grid.
% m order
% L_sig is the band-limit
f_temp=zeros(size(FI_t));
f_temp_neg=zeros(size(FI_t));

if mod(m,2) == 0 % m even
    gm = transpose(P_mat(:,1:m/2))*flm_t; 
    gm_neg = conj(transpose(P_mat(:,1:m/2))*flm_t);
    %assigning gm(theta) to the correct  f(theta, phi)
    for ii=0:m/2-1
       f_temp(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm(ii+1); 
       f_temp_neg(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm_neg(ii+1); 
    end
else
     gm = transpose(P_mat(:,1:(m+1)/2))*flm_t; 
    gm_neg = conj(transpose(P_mat(:,1:(m+1)/2))*flm_t);
    %assigning gm(theta) to the correct  f(theta, phi)
    for ii=0:(m+1)/2-1
       f_temp(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm(ii+1); 
       f_temp_neg(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm_neg(ii+1); 
    end
end




% no need of condition on m here because the spatial elimination is not run
% for m=0. (Verified. 15/11/2013)

f_t = f_t - f_temp.*exp(1i*m*FI_t) - f_temp_neg.*exp(-1i*m*FI_t);

end







function flm_t  = nsht_inversion(P_mat, fft_v,m,L_sig,n_rad,lambda_l,lambda_n)
% nsht_inversion - when passed with fft vector and order m,
% band-limit L, returns the spherical harmonic coefficients. This function
% would perform matrix inversion.

% fft_v should have the same structure as of FI
% flm_t is the vector which contains spectral coefficeints [fm^m, f(m+1)^m
% fL^m];

% Ax = b, first form b from fft_vec

b = zeros(1,L_sig-m); 
for el = m:L_sig-1
      b(el-m+1) = fft_v(m+el^2+1);     
end

%keep non-zero values
if mod(m,2) == 0 % m even
    b = b(1:2:L_sig-m);
    P = P_mat(1:2:end,m/2+1:end);
    i = m:2:(L_sig-1);
else
    b = b(2:2:L_sig-m);
    P = P_mat(2:2:end,(m+1)/2+1:end);
    i =(m+1):2:(L_sig-1);
end


LB = diag(i.^2.*(i+1).^2);
N_mat = diag(repmat(n_rad^2*(n_rad+1)^2,length(i),1));

P = transpose(P);
b = transpose(b);
flm_short = (P'*P + lambda_l*LB+ lambda_n*N_mat)\(P'*b);

flm_t = zeros(1,L_sig-m);
if mod(m,2) == 0 % m even
    for i = 1: length(flm_short)
        flm_t((i*2)-1) = flm_short(i);
    end
else
     for i = 1: length(flm_short)
        flm_t((i*2)) = flm_short(i);
    end

end

flm_t = transpose(flm_t);

end

