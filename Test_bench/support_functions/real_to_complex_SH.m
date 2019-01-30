function [complex_coeffs] =real_to_complex_SH(real_coeffs, Lmax)
%converts real even degree SH to complex SH
% Lmax is the max SH degree
%date last changed: 12/12/17
%author:Alice Bates



%if vector
if (size(real_coeffs,2)==1) || (size(real_coeffs,1)==1) 
    m =addmout(Lmax,1);
    complex_coeffs = zeros((Lmax+1)*(Lmax+2)*0.5,1);

    for  i=1: ((Lmax+1)*(Lmax+2)*0.5)
        if m(i) == 0
            complex_coeffs(i) = real_coeffs(i);
        elseif m(i) < 0
            complex_coeffs(i) = 1/sqrt(2)*real_coeffs(i-2*m(i)) + 1i/sqrt(2)*real_coeffs(i);           
        else
            complex_coeffs(i) = (-1)^m(i)/sqrt(2)*real_coeffs(i) + (-1)^(m(i)+1)*1i/sqrt(2)*real_coeffs(i-2*m(i));  
        end
    end

else
%if matrix of coefficient
    complex_coeffs = zeros(size(real_coeffs));
    for j=1:size(real_coeffs,2)
        complex_coeffs(:,j) = real_to_complex_SH(real_coeffs(:,j), Lmax);
    end
end




 
