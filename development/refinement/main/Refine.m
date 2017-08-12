function [ v_new ] = Refine( v, data, params, iter_max, tol, CTF_flag, filename )
%This is the main function for refinement
%   Input
%           v: reference volume.
%           data: experimental projection images
%           params: data parameters
%                   params.N: number of defocus groups
%                   params.defidx: image's defocus group index
%                   params.max_shifts: maximum shifts in x, y directions
%                   params.c: The experimental CTF functions.
%           iter_max: maximum refinement iteration
%           tol: tolerance for convergence of the refinement.
%	    CTF_flag: if the volume v is CTF corrected, CTF_flag = 1. If v is not CTF corrected, CTF_flag=0.	 
%           filename: file name for storing the refinement results.
%   Output:
%           v_new: final refined volume
%%Zhizhen Zhao 02/2014

v_new=zeros(size(v));
L=size(v, 1);
P=size(data, 3);
ref_k=1000; %generate 1000 reference images

N = params.N;
d = params.defidx;
max_shifts=params.max_shifts;
c=params.c;

r_max=floor(L/2)-max_shifts;
class=zeros(P, 1);
corr=zeros(P, 1);
rot=zeros(P, 1);
shifts_b=zeros(P, 2);
iter=0;
norm_diff_ratio = 1.0;
q = qrand(ref_k);

if ~exist('./results','dir')
    if ~mkdir('./results')
        error('Cannot create directory ./results');
    end
end

while (iter<iter_max && norm_diff_ratio>tol )
    filename=sprintf('%s_iter%d', filename, iter+1); %save the results for each iteration
    
    [ref]=cryo_project(real(v), q); %generate references
    ref = permute(ref, [2, 1, 3]);
    %
    for i=1:N
        id=find(d==i);
        if (iter==0 && CTF_flag==0)
            [corr(id), class(id), rot(id), shifts_b(id, :)] = refinement_Bessel(data(:, :, id), ref, ones(L), max_shifts, r_max );
        else
            [corr(id), class(id), rot(id), shifts_b(id, :)] = refinement_Bessel(data(:, :, id), ref, c(:, :, i), max_shifts, r_max );
        end;
        fprintf('\nFinished defocus group %d.', i);
    end
    
    %%%%Align images
    [data2]=refinement_align(rot, shifts_b, data);
    
    R=zeros(3, 3, ref_k);
    for i=1:ref_k
        R(:, :, i)=q_to_rot(q(:, i));
    end;
    inv_rot=zeros(3, 3, length(d));
    for i=1:length(d)
        inv_rot(:, :, i)=R(:, :, class(i)).';
    end;
    
    [~, id_select] = sort(corr, 'descend');
    id_select=id_select(1:floor(0.8*P)); %select images with high normalized cross correlation value with the references.
    [ v_new ] = recon3d_firm_ctf( data2(:, :, id_select),...
        c, d(id_select),inv_rot(:, :, id_select),[], [], 30, [] );
    v_new=real(v_new);
   
    norm_diff_ratio = norm(v_new(:) - v(:))/norm(v(:));
    
    v = v_new;
    
    save(filename, 'v_new', 'class', 'rot', 'shifts_b', 'corr', 'iter');

    iter=iter+1;
    fprintf('\nIteration number: %d', iter);

end
