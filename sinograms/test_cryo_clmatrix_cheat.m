function test_cryo_clmatrix_cheat
%
% TEST_CRYO_CLMATRIX_CHEAT Test consistency between common lines routines.
%
% Compare the output of cryo_clmarix and cryo_clmatrix_cheat.
%
% Yoel Shkolnisky, November 2017

n=65;
K=10;
precision = 'single';
mrcname= fullfile(aspire_root,'./projections/simulation/maps/cleanrib.mrc');
vol1 = ReadMRC(mrcname);
rots = rand_rots(K);
projs = cryo_project(vol1,rots,n,precision);
projs = permute(projs, [2 1 3]);


% for i = 1:size(rots,3)
%     rots(:,:,i) = rots(:,:,i)';
% end
% 
n_r=65;
n_theta=360;
npf=cryo_pft(projs,n_r,n_theta,'single');

[clmatrix_ref,~,~] = clmatrix_cheat(rots,n_theta);
clmatrix = cryo_clmatrix_gpu(npf, -1, 1,1, 0.1);

% Find mismathces
max_angle=5;
[~,matches]=comparecl(clmatrix,clmatrix_ref,360,max_angle);
log_message('Indices with common lines errors (more than %d degress):',max_angle);
count=0;
for k1=1:K-1
    for k2=k1+1:K
        if matches(k1,k2)==0
            log_message('\t (%d,%d) \t clmatrix(%d,%d) \t reference(%d,%d)',...
                k1,k2,clmatrix(k1,k2),clmatrix(k2,k1),...
                clmatrix_ref(k1,k2),clmatrix_ref(k2,k1));
            count=count+1;
        end
    end
end
log_message('Total errors: %d pairs',count);
