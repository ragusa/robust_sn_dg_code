function [T, S, F, D, M, Sigma, C]=build_matrices(porder,tot,sca,dx,SNQ,gamma_0,delta_0,logi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute elementary matrices
% mt = mass matrix for total xs
% ms = mass matrix for scattering xs
% g  = gradient matrix
% e =  edge matrix e{1} for mu>0, e{2} for mu<0
[mt ,ms ,g ,e ]=compute_elem1(porder,tot,sca,dx,gamma_0,delta_0,logi );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute global matrices
% T= transport matrix acting on angular fluxes
% S= scattering operator acting on angular fluxes
% F= transport matrix acting on scalar fluxes
% D= discrete to moment matrix
% Sigma= scattering matrix acting on moments
% M= moment to discrete matrix
% note: M.Sigma.D = S
% C= matrix acting on angular fluxes and producing the current

[T ,S ,F ,D ,M ,Sigma ,C ]=compute_T1(porder,SNQ,mt,ms,g,e ,dx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

