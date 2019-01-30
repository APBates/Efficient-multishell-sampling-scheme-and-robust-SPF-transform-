function [EM,EL]=addmout(L,evens)
% addmout.m returns spherical harmonic degree and order indexing arrays.
% For arrays where m=-l:l
%
% INPUT:
%
% L        Maximal degree of the expansion (bandwidth)
% evens    0 For all degrees [default]
%          1 Only do the even degrees
%
% OUTPUT:
%
% EM       Vector of orders  m involved in the expansion, [0 -101 -2-1012] 
% EL       Vector of degrees l involved in the expansion, [0  111  2 2222]
%
%
% Date last changed: 14/02/18

matr=(repmat(0:evens+1:L,2,1)'*diag([-1 1]))';
EM=matranges(matr(:)')';
twolp=2*(0:evens+1:L)+1;
EL=gamini(0:evens+1:L,twolp)';



