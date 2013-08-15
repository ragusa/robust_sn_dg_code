function phi=compu_phi(psi,ndof,SNQ);
% compute scalar flux and current

% initialize
phi=zeros(ndof,2);

% compute moments
for idir=1:SNQ.sn
    i1=(idir-1)*ndof + 1;
    i2=(idir  )*ndof   ;
    phi = phi + SNQ.w(idir)* kron( psi(i1:i2),  [ 1 SNQ.mu(idir)] );
end
