function [phi]=direct_solve_phi(F,T,D,qva_sw,qsa);

% direct solve for phi
bb = D*inv(T) * qva_sw;
aa = D*inv(T) * qsa;
n0=length(F);

phi_ = inv(speye(n0)-F) * ( bb + aa);

phi(:,1) = phi_(1:n0/2);
phi(:,2) = phi_(n0/2+1:end);
