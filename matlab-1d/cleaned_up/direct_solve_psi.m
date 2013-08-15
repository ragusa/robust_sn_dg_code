function [phi,psi_]=direct_solve_psi(T,S,C,q,ndof,SNQ);

% direct solve for the angular flux
psi = inv(T-S+C) * ( q );

% postprocess to obtain the scalar flux and the current (first moment of
% psi)
phi=compu_phi(psi,ndof,SNQ);

if(nargout==2)
    psi_=psi;
end
