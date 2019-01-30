function[Elmn] = fSPFt(E,L,N,x,zeta,lambda_l, lambda_n)
%fSPFt performs the forward transform in the SPF basis to calculate the SPF
%coefficients. This version allows angular and radial regularisation 
% Inputs:
% E          the diffusion signal stored as a N x (Lmax+1)*Lmax*0.5 matrix.
% i.e. each row contains samples in one shell, rows padded with zeros at
% end for shells with less samples
% L          a vector containing the SH band-limit at different shells.
% This is the max SH degree.
% N          the radial truncation order = number of shells. Max radial
% order = N-1
% x         Gauss-Laguerre quadrature roots 
% zeta       scale factor
% Output:
% Elmn      the coefficients in the SPF basis stored as a N x Lmax^2
% matrix
%author: Alice Bates
%date: 12/02/17

Lmax = max(L);
Elmn_3D = zeros(N+1,(Lmax+1)^2,N+1);

%going through each shell and perform SHT
for q=1:N+1
    Elm = zeros(N+1,(Lmax+1)^2);
    E_shell = E(q,1:(L(q)+1)*(L(q)+2)*0.5);
    w =0.5*zeta^(1.5)*(gamma(N+2.5)*x(q)*exp(x(q)))/( factorial(N+1)*(N+2)^2*laguerreL(N+2,0.5,x(q))^2 );
    n_vec = 0:N;
    
    %perform SHT for every radial degree n for the current shell q
    Elm(:,1:(L(q)+1)^2)=nsht_forward_reg_all_n(E_shell,L(q),q,lambda_l,lambda_n,n_vec);

    
    %work out value of radial function at roots
    constant_n = ( (2* factorial(n_vec))./(zeta^1.5*gamma(n_vec+1.5)) ).^0.5;
    R =  transpose(constant_n.*exp(-x(q)./2).*laguerreL(n_vec,0.5,x(q)));  
    
    Elmn_3D(:,:,q) = repmat(w*R,1,(Lmax+1)^2).*Elm;
end

%sum over each shell
Elmn = sum(Elmn_3D,3);