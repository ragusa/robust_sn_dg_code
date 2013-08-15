function [mtot,msca,gr,e]=compute_elem1(porder,tot,sca,dx,gamma_0,delta_0,logi_std_up)

% # of ncells
ncells = length(dx);
% ndofs
ndofs=(porder+1)*ncells;

% enforce standard upwind scheme
if(logi_std_up)
    fu=1; fd=0; delta=0;
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
            [gamma,delta]=comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0);
            fu=(1+gamma)/2;
            fd=(1-gamma)/2;

            % interior edges
            if(iel==1)
                e{1}(istart+1,istart +1) =  fu;
                e{1}(istart+1,istart2  ) =  fd;
            end
            if(iel>1)&&(iel<ncells)
                e{1}(istart  ,istart0+1) = -fu;
                e{1}(istart  ,istart   ) = -fd;
                e{1}(istart+1,istart +1) =  fu;
                e{1}(istart+1,istart2  ) =  fd;
            end
            if(iel==ncells)
                e{1}(istart  ,istart0+1) = -fu;
                e{1}(istart  ,istart   ) = -fd;
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
            [gamma,delta]=comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0,ncells);
            fu=(1+gamma)/2;
            fd=(1-gamma)/2;

            % interior edges
            if(iel==1)
                e{2}(istart+1,istart +1) =  fd;
                e{2}(istart+1,istart2  ) =  fu;
            end
            if(iel>1)&&(iel<ncells)
                e{2}(istart  ,istart0+1) = -fd;
                e{2}(istart  ,istart   ) = -fu;
                e{2}(istart+1,istart +1) =  fd;
                e{2}(istart+1,istart2  ) =  fu;
            end
            if(iel==ncells)
                e{2}(istart  ,istart0+1) = -fd;
                e{2}(istart  ,istart   ) = -fu;
            end
            % outgoing face
            if(iel==1), e{2}(istart,istart ) = -1; end
        end

        % current edge matrix
        % edge matrix,
        if(iel~=1)
            e{3}(istart,istart0+1 ) = -1 *delta;
            e{3}(istart,istart    ) =  1 *delta;
        end
        if(iel~=ncells)
            e{3}(istart+1,istart+1 ) =  1 *delta;
            e{3}(istart+1,istart2  ) = -1 *delta;
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
            [gamma,delta]=comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0);
            fu=(1+gamma)/2;
            fd=(1-gamma)/2;

            % interior edges
            if(iel==1)
                e{1}(iel,iel  ) =   fu;
                e{1}(iel,iel+1) =   fd;
            end
            if(iel>1)&&(iel<ncells)
                e{1}(iel,iel-1) =  -fu;
                e{1}(iel,iel  ) =  -fd+fu;
                e{1}(iel,iel+1) =      fd;
            end
            if(iel==ncells)
                e{1}(iel,iel-1) =  -fu;
                e{1}(iel,iel  ) =  -fd;
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
            [gamma,delta]=comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0,ncells);
            fu=(1+gamma)/2;
            fd=(1-gamma)/2;

            % interior edges
            if(iel==1)
                e{2}(iel,iel  ) =   fd;
                e{2}(iel,iel+1) =   fu;
            end
            if(iel>1)&&(iel<ncells)
                e{2}(iel,iel-1) =  -fd;
                e{2}(iel,iel  ) =  -fu+fd;
                e{2}(iel,iel+1) =      fu;
            end
            if(iel==ncells)
                e{2}(iel,iel-1) =  -fd;
                e{2}(iel,iel  ) =  -fu;
            end
            % outgoing face
            if(iel==1), e{2}(iel,iel) = e{2}(iel,iel)-1; end
         end

        % current edge matrix
        switch iel
            case 1
                e{3}(iel,iel:iel+1) = [1 -1] *delta;
            case ncells
                e{3}(iel,iel-1:iel) = [-1 1] *delta;
            otherwise
                e{3}(iel,iel-1:iel+1) = [-1 2 -1] *delta;
        end
    end

end


