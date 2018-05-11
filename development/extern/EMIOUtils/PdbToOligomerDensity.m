% PdbToOligomerDensity
% Read a pdb file and construct a 3D density map that is Gaussian filtered
% and downsampled to the given resolution.

name='3NAF.pdb';
fc=1/17; % Inverse A for Gauss filter
res=2;  % must be an integer
border=1.2; % fraction of oversampling

disp('Reading .pdb file');
% ReadPDB;  % returns coords, rots1
[coords rots1]=ReadPDB(name);

nrots=size(rots1,3);  % number of rotation matrices
if nrots==0
    rots1=eye(3,3);
    rots1(:,4)=zeros(3,1);
    nrots=1;
end;
%%
nj=size(coords,2);
cors=[];
for i=1:nrots
    cors=[cors round(rots1(:,1:3,i)*coords+repmat(rots1(:,4,i),1,nj))];
end;
maxs=max(cors')';
mins=min(cors')';
mids=(maxs+mins)/2;

n1=ceil(max(maxs-mins)+1);  % size of minimum bounding cube
n=NextNiceNumber(n1*1.1,5,res*2); % must be a multiple of res*2 to allow downsampling.
n
vol=zeros(n,n,n);
ctr=n/2+1-mids;
%%
cors=cors+repmat(ctr,1,nj*nrots);

%%
disp('Inserting volume');
for j=1:size(cors,2)
    vol(cors(1,j),cors(2,j),cors(3,j))...
    =1+vol(cors(1,j),cors(2,j),cors(3,j));
end;

%%
disp('filtering...')
fc=1/17;  % 17 A filter
fv=SharpFilt(vol,fc,fc/10);
fvd=Downsample(fv,n/res);
ShowSections(fvd);


WriteMRC(fvd,res,'3NAF.mrc');

