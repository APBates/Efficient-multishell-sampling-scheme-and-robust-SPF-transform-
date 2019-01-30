function y = mBessel_ratio(n,x)
% y = mBessel{n}(x)/mBessel{n-1}(x) = besseli(n,x)./besseli(n-1,x)
% Fast evaluation using the Perron's Continued Fraction equation.
% For more details see: http://www.ams.org/journals/mcom/1978-32-143/S0025-5718-1978-0470267-9/

y = x./( (2*n + x) - ( 2*x.*(n+1/2)./ ( 2*n + 1 + 2*x - ( 2*x.*(n+3/2)./ ( 2*n + 2 + 2*x - ( 2*x.*(n+5/2)./ ( 2*n + 3 + 2*x ) ) ) ) ) ) );
end 