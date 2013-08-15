function [mtot,msca,gr,e]=compute_elem1(porder,tot,sca,dx,gamma_0,delta_0,logi_std_up)

% # of ncells
ncells = length(dx);
% ndofs
ndofs=(porder+1)*ncells;

% enforce standard upwind scheme
if(logi_std_up)
    fu=1; fd=0; deltaL=0; deltaR=0;
end

if(porder==1)
    % mass matrix
    m = [2 1; 1 2]/3;
    % gradient matrix
    g = [1 1 ; -1 -1]/2;
else
    m=[2];
    g=[0];
end
% jacobian of the affine transformation 
jac  = dx(1:ncells)/2;

% global mass matrices
mtot =  kron( sparse( diag(tot(1:ncells).*jac) ), sparse(m) );
for ani=1:length(sca)
    msca{ani} =  kron( sparse( diag(sca{ani}(1:ncells).*jac) ), sparse(m) );
end
% global gradient matrix
gr = kron( speye(ncells), sparse(g) );

% edge matrix
e{1} = spalloc(ndofs,ndofs,ndofs);
e{2} = spalloc(ndofs,ndofs,ndofs);
% current edge matrix
e{3} = spalloc(ndofs,ndofs,ndofs);

if(porder==1)
    % build edge matrix
    for iel=1:ncells
        istart = (iel-1)*2+1;
        iel0 = iel-1;
        istart0 = (iel0-1)*2+1;
        iel2 = iel+1;
        istart2 = (iel2-1)*2+1;

        % edge matrix, mu>0
        % ~~~~~~~~~~~~~~~~~
        if(logi_std_up)
            if(iel~=1),      e{1}(istart  ,istart0+1) = -fu; end
            e{1}(istart+1,istart+1 ) =  fu;
        else
            [gammaL,gammaR,deltaL,deltaR]=comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0,ncells);
            fuL=(1+gammaL)/2;
            fdL=(1-gammaL)/2;
            fuR=(1+gammaR)/2;
            fdR=(1-gammaR)/2;

            % interior edges
            if(iel==1)
                e{1}(istart+1,istart +1) =  fuR;
                e{1}(istart+1,istart2  ) =  fdR;
            end
            if(iel>1)&(iel<ncells)
                e{1}(istart  ,istart0+1) = -fuL;
                e{1}(istart  ,istart   ) = -fdL;
                e{1}(istart+1,istart +1) =  fuR;
                e{1}(istart+1,istart2  ) =  fdR;
            end
            if(iel==ncells)
                e{1}(istart  ,istart0+1) = -fuL;
                e{1}(istart  ,istart   ) = -fdL;
            end
            % outgoing face
            if(iel==ncells), e{1}(istart+1,istart+1 ) = 1; end
        end

        % edge matrix, mu<0
        % ~~~~~~~~~~~~~~~~~
        if(logi_std_up)
            e{2}(istart  ,istart   ) = -fu;
            if(iel~=ncells), e{2}(istart+1,istart2  ) =  fu; end
        else
            [gammaL,gammaR,deltaL,deltaR]=comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0,ncells);
            fuL=(1+gammaL)/2;
            fdL=(1-gammaL)/2;
            fuR=(1+gammaR)/2;
            fdR=(1-gammaR)/2;

            % interior edges
            if(iel==1)
                e{2}(istart+1,istart +1) =  fdR;
                e{2}(istart+1,istart2  ) =  fuR;
            end
            if(iel>1)&(iel<ncells)
                e{2}(istart  ,istart0+1) = -fdL;
                e{2}(istart  ,istart   ) = -fuL;
                e{2}(istart+1,istart +1) =  fdR;
                e{2}(istart+1,istart2  ) =  fuR;
            end
            if(iel==ncells)
                e{2}(istart  ,istart0+1) = -fdL;
                e{2}(istart  ,istart   ) = -fuL;
            end
            % outgoing face
            if(iel==1), e{2}(istart,istart ) = -1; end
        end

        % current edge matrix
        % edge matrix,
        if(iel~=1)
            e{3}(istart,istart0+1 ) = -1 *deltaL;
            e{3}(istart,istart    ) =  1 *deltaL;
        end
        if(iel~=ncells)
            e{3}(istart+1,istart+1 ) =  1 *deltaR;
            e{3}(istart+1,istart2  ) = -1 *deltaR;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
else % (porder =0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

    % build edge matrix
    for iel=1:ncells

        iel0 = iel-1;
        iel2 = iel+1;

        % edge matrix, mu>0
        % ~~~~~~~~~~~~~~~~~
        if(logi_std_up)
            if(iel~=1),      e{1}(iel,iel-1) = -fu; end
            e{1}(iel,iel) =  fu;
        else
            [gammaL,gammaR,deltaL,deltaR]=comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0,ncells);
            fuL=(1+gammaL)/2;
            fdL=(1-gammaL)/2;
            fuR=(1+gammaR)/2;
            fdR=(1-gammaR)/2;

            % interior edges
            if(iel==1)
                e{1}(iel,iel  ) =   fuR;
                e{1}(iel,iel+1) =   fdR;
            end
            if(iel>1)&(iel<ncells)
                e{1}(iel,iel-1) =  -fuL;
                e{1}(iel,iel  ) =  -fdL+fuR;
                e{1}(iel,iel+1) =       fdR;
            end
            if(iel==ncells)
                e{1}(iel,iel-1) =  -fuL;
                e{1}(iel,iel  ) =  -fdL;
            end
            % outgoing face
            if(iel==ncells), e{1}(iel,iel) = e{1}(iel,iel)+1; end
        end

        % edge matrix, mu<0
        % ~~~~~~~~~~~~~~~~~
         if(logi_std_up)
            if(iel~=ncells),      e{2}(iel,iel+1) = fu; end
            e{2}(iel,iel) =  -fu;
        else
            [gammaL,gammaR,deltaL,deltaR]=comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0,ncells);
            fuL=(1+gammaL)/2;
            fdL=(1-gammaL)/2;
            fuR=(1+gammaR)/2;
            fdR=(1-gammaR)/2;

            % interior edges
            if(iel==1)
                e{2}(iel,iel  ) =   fdR;
                e{2}(iel,iel+1) =   fuR;
            end
            if(iel>1)&(iel<ncells)
                e{2}(iel,iel-1) =  -fdL;
                e{2}(iel,iel  ) =  -fuL+fdR;
                e{2}(iel,iel+1) =       fuR;
            end
            if(iel==ncells)
                e{2}(iel,iel-1) =  -fdL;
                e{2}(iel,iel  ) =  -fuL;
            end
            % outgoing face
            if(iel==1), e{2}(iel,iel) = e{2}(iel,iel)-1; end
         end

        % current edge matrix
        switch iel
            case 1
                e{3}(iel,iel:iel+1) = [1 -1] *deltaR;
            case ncells
                e{3}(iel,iel-1:iel) = [-1 1] *deltaL;
            otherwise
                e{3}(iel,iel-1:iel+1) = [-1 1 0] *deltaL + [0 1 -1] *deltaR; % [-1 2 -1] *delta;
        end
    end

end


