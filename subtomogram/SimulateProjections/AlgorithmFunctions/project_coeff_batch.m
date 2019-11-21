function [coeff_pos_k] = project_coeff_batch(CR, Crm, InnerP, L, N, ns, nn)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% fprintf('project_coeff_batch\n'); fprintf(datestr(now));
Psilm = cell(1,L); % index l*m, inside cell rotation*q
batch_size = N;%floor(1000/ns);
batch_num = 0;%floor(N/batch_size);
batch_rem = N;%mod(N, batch_size);
nls = zeros(L,1);
for ll = 1:L
    nls(ll) = size(InnerP{ll,ll},2);
    Psilm{1, ll} = zeros(N*ns, nls(ll));
end

for ll = 0:L-1
    for i = 1:batch_num
        temp = permute(CR{ll+1,1}(:,:,(i-1)*batch_size+1:i*batch_size),[3 1 2]);
        temp = mat2cell(reshape(temp,batch_size,[]),batch_size,size(temp,2)*ones(1,ll+1));
        temp = cellfun(@(x,y,z) bsxfun(@times,x*y,z), temp, InnerP(ll+1,1:ll+1), Crm(ll+1,1:ll+1),'UniformOutput',false);
        temp = cellfun(@(x) reshape(permute(x,[3 1 2]),batch_size*ns,[]), temp,'UniformOutput',false);
        Psilm(1,1:ll+1) = cellfun(@(x,y) x+y, ...
            Psilm(1,1:ll+1), temp, 'UniformOutput', false);
    end
    
    if batch_rem ~= 0
        temp = permute(CR{ll+1,1}(:,:,end-batch_rem+1:end),[3 1 2]);
        temp = mat2cell(reshape(temp,batch_rem,[]),batch_rem,size(temp,2)*ones(1,ll+1));
        temp = cellfun(@(x,y,z) bsxfun(@times,x*y,z), temp, InnerP(ll+1,1:ll+1), Crm(ll+1,1:ll+1),'UniformOutput',false);
        temp = cellfun(@(x) reshape(permute(x,[3 1 2]),batch_rem*ns,[]), temp,'UniformOutput',false);
        Psilm(1,1:ll+1) = cellfun(@(x,y) x+y, ...
            Psilm(1,1:ll+1), temp, 'UniformOutput', false);
    end
end

coeff_pos_k = zeros(nn,N*ns);
iter = 0;
for ll = 1:L
%     nl = size(Psilm{end,ll},2);
%     temp = reshape(cell2mat(Psilm(ll:end,ll)),N*ns,L-ll+1,nl);
%     temp = squeeze(sum(temp,2));
    nl = nls(ll);
    coeff_pos_k(iter+1:iter+nl,:) = transpose(Psilm{1,ll});
    iter = iter+nl;
end

end

