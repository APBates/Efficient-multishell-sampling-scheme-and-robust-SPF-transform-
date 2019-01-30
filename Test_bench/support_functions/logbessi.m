function out = logbessi(n,x)
% Computes log(I_n(input)), where I_n is the n-order modified Bessel
% function of the first kind

% Justin Haldar 07/26/2012
out = abs(real(x))+log(besseli(n,x,1));

return;