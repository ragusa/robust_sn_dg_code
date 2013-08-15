function [L,S,F,D,SS,myms,C]=compute_T1(porder,snq,mt,ms,gg,e,dx)

% compute the sum of the SN weights
sw = sum(snq.w);

% dimension of the linear system
ncells = length(dx);
n = ncells*(porder+1);
L = spalloc(n,n,snq.sn*ncells*(porder+1)^2);
% twice as big for phi and J
F = zeros(2*n,2*n);
% anisotropic scattering matrix
z=ms{1}*0;
if(length(ms)==2)
    myms = [ ms{1} z; z ms{2}];
else
    myms = [ ms{1} z; z z ];
end

W  = kron(sparse(diag(snq.w)),speye(n));            % angular weighting
SS = kron( [ ones(snq.sn,1) , snq.mu'], speye(n) ); % spherical harmonics
D  = sparse(SS'*W);
SS=SS/sw;

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WHAT IS BELOW IS FOR THE DIRECT SOLVER BASED ON THE FLUX MOMENTS, NOT
    % FOR THE DIRECT SOLVER BASED ON THE ANGULAR FLUXES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the other way is to perform the sweep and collapse in angle to
    % obtain directly an error matrix in the scalar errors
    Md=[ speye(n)    snq.mu(idir)*speye(n)];
    %     Dd=[ snq.w(idir)*speye(n) ; snq.w(idir)*snq.mu(idir)*speye(n)];
    Dd=snq.w(idir) * Md';
    %     F = F + Dd * inv(Ld)* Md * (myms/sw);
    %     F = F + Dd * inv(Ld)* ( Md * myms - [z e{3}] ) / sw;
    mycur = [z z;z e{3}];
    F = F + Dd * inv(Ld)* Md * ( myms - mycur ) / sw;
end




