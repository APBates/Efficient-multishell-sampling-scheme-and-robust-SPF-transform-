function [real_coeffs] = complex_to_real_SH(complex_coeffs, Lmax)
%converts real even degree SH to complex SH
%Lmax is the max SH degree
%date last changed: 6/11/17
%author:Alice Bates

%if row or column vector
if (size(complex_coeffs,2)==1) || (size(complex_coeffs,1)==1) 
    m =addmout(Lmax,1);
    real_coeffs = zeros((Lmax+1)*(Lmax+2)*0.5,1);

    for  i=1: ((Lmax+1)*(Lmax+2)*0.5)
        if m(i) == 0
            real_coeffs(i) = complex_coeffs(i);
        elseif m(i) < 0
            real_coeffs(i) = sqrt(2)*imag(complex_coeffs(i));         
        else
            real_coeffs(i) = sqrt(2)*real(complex_coeffs(i-2*m(i)));
        end
    end
    real_coeffs = real(real_coeffs);
else
% if matrix of coefficient
    real_coeffs = zeros(size(complex_coeffs));
    for j=1:size(complex_coeffs,2)
        real_coeffs(:,j) = complex_to_real_SH(complex_coeffs(:,j), Lmax);
    end
end


