function [ODF_lm] = ODF_SH_coefficients(zeta, Elmn,Lmax,N)
%calculates the ODF from the SPF coefficients - may need to convert the SPF
%coefficients to real coefficients
%author: Alice Bates
%date last changed: 23/02/18

l_vec = [0,repmat(2,1,5),repmat(4,1,9),repmat(6,1,13),repmat(8,1,17),repmat(10,1,21)];

ODF_lm = zeros((Lmax+1)*(Lmax+2)*0.5,1);
%for each coefficient

constant_n = zeros(1,N);
for ind=1:N
    constant_n(ind) = ( (2* factorial(ind))/(zeta^1.5*gamma((ind)+1.5)) )^0.5;
end



for index=1:(Lmax+1)*(Lmax+2)*0.5
    l = l_vec(index);
    %m = index-1-0.5*l^2-0.5*l;
    if l==0
        ODF_lm(index) = 1/sqrt(4*pi);
    else
        sum_n=0;
        for n=1:N
            sum_i = 0; 
            for i=1:n
                 a = nchoosek(n+0.5,n-i);
                 if isempty(a)
                    a=0;
                end
                sum_i=sum_i+(-1)^i*a*2^i/i;
            end
           sum_n = sum_n + constant_n(n)*Elmn(n+1,index)*sum_i;               
        end
         ODF_lm(index) = -1/(8*pi)*legendreP(l,0)*(-l)*(l+1)*sum_n;
    end
end

