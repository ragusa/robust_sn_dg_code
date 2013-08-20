function [R]=build_residual(porder,tot,sca,dx,snq,gamma_0,delta_0,logi_std_up,psi,q)

moments = compu_phi(psi,ndof,snq);
phi = moments(:,1);

[mt ,ms ,gg ,e ]=compute_elem1(porder,tot,sca,dx,gamma_0,delta_0,logi_std_up,phi );

% compute the sum of the SN weights
sw = sum(snq.w);

% dimension of the linear system
ncells = length(dx);
n = ncells*(porder+1);
L = spalloc(n,n,snq.sn*ncells*(porder+1)^2);
% anisotropic scattering matrix
z=ms{1}*0;
if(length(ms)==2)
    myms = [ ms{1} z; z ms{2}];
else
    myms = [ ms{1} z; z z ];
end

% W  = kron(sparse(diag(snq.w)),speye(n));            % angular weighting
% SS = kron( [ ones(snq.sn,1)/sw , snq.mu'], speye(n) ); % spherical harmonics
% SS=SS/sw;

% build the action of all transport sweeps on the scattering term
ndir = snq.sn;
for idir=1:ndir
    m=snq.mu(idir); % shortcut
    % choose the edge contribution based on the sweeping order (L to R or R to L)
    if(m>0)
        edg=e{1};
    else
        edg=e{2};
    end

    i1 = (idir-1)*n+1;
    i2 = idir*n;

    % transport matrix for this direction
    % = mu*gradient + total mass + edge
    % note that the jacobian dx/2 is already included in the mass matrices
    % and that it does not appear in the gradient and edge matrices for 1D
    Ld = ( m*gg + mt ) + m*edg;

    % one way is to build the entire SN matrix T
    % and the scattering matrix S
    L(i1:i2,i1:i2) = Ld;
    S(i1:i2,:) = kron(snq.w, (ms{1}) ) /sw;
    if(length(ms)==2)
        S(i1:i2,:) = S(i1:i2,:) + kron(...
            snq.mu .* snq.w, snq.mu(idir)*ms{2} ) /sw;
    end

    % current term
    %    C(i1:i2,:) =  kron( snq.mu .* snq.w,  snq.w(idir)*snq.mu(idir)*e{3} ) /sw;
    C(i1:i2,:) =  kron( snq.mu .* snq.w,  snq.mu(idir)*e{3} ) /sw;
    %     C(i1:i2,:) =  kron( snq.mu .* snq.w,  e{3} ) /sw;

end

% direct solve for the angular flux
(T-S+C)*psi - q; %(qva/sw+qsa) ;

