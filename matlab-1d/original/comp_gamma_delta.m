function [gamma,delta] = comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0,n)

% take sigma_s from the element where the basis function is evaluated
sig = sca{1}(iel);

%%%
%%% uncomment below for harmonic average
%%% 
% % if(nargin==6)
% %     dir=+1;
% % else
% %     dir=-1;
% % end
% % % previous and next element number
% % iel0=iel-1;
% % iel2=iel+1;
% % if(dir>0)
% %     i=iel0;
% %     if(i==0), i=iel;end % std boundary terms
% % else
% %     i=iel2;
% %     if(i==n+1), i=iel;end % std boundary terms
% % end
% % i_neigh=i;
% % i=iel;
% % isig_i     =  1/sca{1}(i) ;
% % isig_neigh =  1/sca{1}(i_neigh) ;
% % sig = 2/( isig_i + isig_neigh );

% comute gamma
gamma = gamma_0 / ( max( gamma_0, sig*sum(dx) ));
gamma = min(1,gamma);

% compute delta
delta = delta_0/gamma*(1-gamma);
