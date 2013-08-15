clear all; % clears all variables from memory
close all; % closes all plotting figures
clc;       % clear console screen

% select spatial approx order (0 or 1)
porder=0;
% select angular approx (must eb an even number)
sn=8;

% sanity checks
if porder~=0 || porder~=1
    error('porder must be 0 or 1')
end
if mod(sn,2)~=0
    error('sn must be even')
end

% load the angular quadrature
[SNQ] = loadquadrature(sn);
% sum of the weights
sw = sum(SNQ.w);

% load data
dataID=5;
[tot,sca,qva,qsa,dx] = loadmydata(dataID,porder,SNQ);
ndof=length(dx)*(porder+1);
% constants gamma_0 and delta_0, as in the paper
delta_0 = 1;
gamma_0 = 2.;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T= transport matrix acting on angular fluxes
% S= scattering operator acting on angular fluxes
% F= transport matrix acting on scalar fluxes
% D= discrete to moment matrix
% Sigma= scattering matrix acting on moments
% M= moment to discrete matri
% note: M.Sigma.D = S
% C= matrix acting on angular fluxes and producing the current

% build standard upwind matrices
% note: Cup=0 in this case
% logical is set to true, so gamma=1 and delta=0 are imposed in the routine
[Tup, S, Fup, D, M, Sigma, Cup] = build_matrices(porder,tot,sca,dx,SNQ,1.,0.,true);

% direct solve for psi: standard upwind method
tic
phiD_up = direct_solve_psi(Tup,S,Cup,qva/sw+qsa,ndof,SNQ);
toc

% only pass the first moment, i.e., the scalar flux
figID=1; figure(figID); %clf;
myplot2(figID,phiD_up(:,1),porder,dx,'b+-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
porder=0;
[tot,sca,qva,qsa,dx] = loadmydata(dataID,porder,SNQ);
ndof=length(dx)*(porder+1);

% build new scheme matrices
[T,S_,F,D_,M_,Sigma_,C]=build_matrices(porder,tot,sca,dx,SNQ,gamma_0,delta_0,false);
clear S_ D_ M_ Sigma_ ; % they are the same as S, D, M, Sigma 

% direct solve for psi: new method
tic
[phiD_new,psi] = direct_solve_psi(T,S,C,qva/sw+qsa,ndof,SNQ);
toc

%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%

myplot2(figID,phiD_new(:,1),porder,dx,'c+-');

% legend('Reduced Upwind','Standard Upwind','Location','Best');
xlabel('position','FontSize',12);
ylabel('Scalar flux','FontSize',12);

disp('Done !!!!');
