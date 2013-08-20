function [gamma_left,gamma_right,delta_left,delta_right] = comp_gamma_delta(sca,tot,dx,iel,gamma_0,delta_0,n)

im1 = max(iel-1,1);
ip1 = min(iel+1,n);

sca_i   = sca{1}(iel);
sca_im1 = sca{1}(im1);
sca_ip1 = sca{1}(ip1);

tot_i   = tot(iel);
tot_im1 = tot(im1);
tot_ip1 = tot(ip1);

g0=0.1;g1=g0;

% compute the left-face gamma
if(sca_i>1e-10)
    gamma_1 = g0/max(g1, max(sqrt( tot_i  /sca_i   - 1 ), sca_i  *dx(iel) ));
else
    gamma_1 = 1;
end
if(sca_im1>1e-10)
    gamma_2 = g0/max(g1, max(sqrt( tot_im1/sca_im1 - 1 ), sca_im1*dx(im1) ));
else
    gamma_2 = 1;
end
% [gamma_1 gamma_2]
gamma_left = max(gamma_1,gamma_2);
gamma_left = 2/(1/gamma_1+1/gamma_2);

% compute the right-face gamma
if(sca_ip1>1e-10)
    gamma_2 = g0/max(g1, max(sqrt( tot_ip1/sca_ip1 - 1 ), sca_ip1*dx(ip1) ));
else
    gamma_2 = 1;
end
% [gamma_1 gamma_2]
gamma_right = max(gamma_1,gamma_2);
gamma_right = 2/(1/gamma_1+1/gamma_2);

% % if(iel==n)
% %     gamma_right
% %     gamma_right=1;
% % end
% % 
% % if(iel<n/2)
% %     gamma_left =1;
% %     gamma_right=1;
% % else
% %     gamma_left =0.015/100;
% %     gamma_right=0.015/100;
% % end

% compute delta's
if(gamma_left>1e-10)
    delta_left  = delta_0/gamma_left *(1-gamma_left );
else
    delta_left=1e8;
end
if(gamma_right>1e-10)
    delta_right = delta_0/gamma_right*(1-gamma_right);
else
    delta_right=1e8;
end

% fid=fopen('gamma_new.txt','a');
% fprintf(fid,'%g %g ',gamma_left,gamma_right);

% old version
% take sigma_s from the element where the basis function is evaluated
sig = sca{1}(iel);
gamma  = gamma_0 / ( max( gamma_0, sig      *sum(dx) ));
delta = delta_0/gamma*(1-gamma);

gamma_left=gamma ; gamma_right=gamma;
delta_left=delta ; delta_right=delta;

% fprintf(fid,' %g \n',gamma);
% fclose(fid);


% tot
% [tot_im1 tot_i sca_ip1]
% sca{1}
% [sca_im1 sca_i sca_ip1]
% 
% [gamma_left,gamma_right,delta_left,delta_right]
% gamma_left=1;
% gamma_right=1;
% delta_left=0;
% delta_right=0;


% % %%%
% % %%% uncomment below for harmonic average
% % % %%% 
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
% % 
% % sig = max( sca{1}(i), sca{1}(i_neigh) ); 
% % % compute gamma
% % gamma = gamma_0 / ( max( gamma_0, sig*sum(dx) ));
% % %% gamma = min(1,gamma);
% % 
% % % compute delta
% % delta = delta_0/gamma*(1-gamma);
