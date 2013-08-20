function [tot,sca,qva,qsa,dx,ndof]=loadmydata(dataID,porder,SNQ)

% sanity checks
if porder~=0 && porder~=1
    error('porder must be 0 or 1')
end

sn=SNQ.sn;
mu=SNQ.mu;

switch dataID
    case 1 % problem A from paper
        % number of cells (equal width for each cell)
        reg_cell = [ 40 ];
        % sigt/sigs per region
        sigt=[100];
        sigs=[100];
        % volumetric source value, per region
        qv=[0.01];
        % domain total length
        L =10;
        % incoming flux values
        inc(1:sn)   = 0;

    case 2 % problem B from paper
        % number of cells (equal width for each cell)
        reg_cell = [ 10 ];
        % sigt/sigs per region
        epsilon=1e-4;
        sigt=[100];
        sigs=[100];
        % volumetric source value, per region
        qv=[0];
        % domain total length
        L =10;
        % incoming flux values
        inc(1:sn)   = 0;
        inc(sn/2+1:end) = 1;

    case 3 % problem C from paper
        % number of cells (equal width for each cell)
        reg_cell = [ 50 ];
        % sigt/sigs per region
        epsilon=1e-4;
        sigt=[100];
        sigs=[100];
        % volumetric source value, per region
        qv=[0];
        % domain total length
        L =10;
        % incoming flux values
        inc(1:sn)   = 0;
        inc(sn)=1/SNQ.w(sn);

    case 4 % problem D from paper
        % number of cells (equal width for each cell)
        reg_cell = [ 10 ];
        % sigt/sigs per region
        epsilon=1e-4;
        sigt=[1/epsilon];
        sigs=[1/epsilon-epsilon];
        % volumetric source value, per region
        qv=[epsilon];
        % domain total length
        L =1;
        % incoming flux values
        inc(1:sn)   = 0;

    case 5 % 2-region problem
        % number of cells (equal width for each cell)
        reg_cell = [ 50 50 ];
        % sigt/sigs per region
        sigt=[100 1e4];
        sigs=sigt;
        % volumetric source value, per region
        qv=[1e-2 1e-2];
        % domain total length
        L =1;
        % incoming flux values
        inc(1:sn)   = 0;

    case 6 % 3-region problem
        % number of cells (equal width for each cell)
        reg_cell = [ 50 50 10];
        % sigt/sigs per region
        sigt=[100 100 100];
        sigs=[50 50 50]*0;
        % volumetric source value, per region
        qv=[1 2 1];
        % domain total length
        L =1;
        % incoming flux values
        inc(1:sn)   = 0;

    otherwise
        error('case ID unknown');

end

% a simple check
aux = sigs./sigt;
i=find(aux>1);
if(~isempty(i))
    error('sigs cannot be > than sigt \nA problem occured with regions %i',i);
end
% total number of elements
ncells = sum(reg_cell);
dx     = L/ncells*ones(1,ncells);
% initialize to zero the tot and scat XS
tot=[];
sca{1}=[]; % isotropic component of the scattering
for i=1:length(reg_cell)
    tot    = [tot    sigt(i)*ones(1,reg_cell(i))];
    sca{1} = [sca{1} sigs(i)*ones(1,reg_cell(i))];
end
sca{2} = 0.*tot; % set linearly anisotropic component to 0


% FEM representation of the angular source
% angular volumetric source
ss=1.;if(porder==0),ss=2;end
qva=[];
n1=1;
for i=1:length(reg_cell)
    n2=reg_cell(i);
    qva=[qva qv(i)*kron( dx(n1:n2)/2, ss*ones(1,(porder+1)) )];
end
qva= kron( ones(1,sn), qva );
qva = qva';
% angular surfacic source
qsa = 0*qva;

% add the mu vales to the surfacic source
ndof=(porder+1)*ncells;
for k=sn/2+1:sn % positive dir
    qsa(ndof*(k-1)+1)=inc(k);
end
for k=1:sn/2 % negative dir
    qsa(ndof*k)=inc(k);
end
z=kron(abs(mu)',ones(ndof,1));
qsa=qsa.*z;

% ndof
ndof=length(dx)*(porder+1);