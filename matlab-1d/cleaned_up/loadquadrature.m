function [snq]=loadquadrature(sn)

% get gauss-legendre abscissae and weights
[a,b]=GLNodeWt(sn);
% put data in snq struct
snq.sn=sn;
snq.mu=a';
snq.w =b';
