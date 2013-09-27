function [m pxy pfile]=ReadHDFStack(filename)
% Read a particle stack as created by EMAN2.  If desired, return the
% original image coordinates and filenames too.

%  filename='z07_ptcls.hdf';

q=hdf5info(filename);
qi=q.GroupHierarchy.Groups.Groups.Groups;
ni=numel(qi);

% Pick up images
xy=qi(1).Datasets.Dims;
m=zeros([xy ni]);
for i=1:ni
    pathname=['/MDF/images/' num2str(i-1) '/image'];
    m(:,:,i)=hdf5read(filename,pathname);
end;

% Pick up particle file and source coordinates
pxy=zeros(ni,2);
pfile=cell(10,1);
for i=1:ni
    qa=qi(i).Attributes;
    for j=1:numel(qa)
        if strcmp(qa(j).Shortname,'EMAN.ptcl_source_coord')
            pxy(i,:)=qa(j).Value';
        elseif strcmp(qa(j).Shortname,'EMAN.ptcl_source_image');
            pfile{i}=qa(j).Value.Data;
        end;
    end;
end;

% % give all the attributes for the first image.
% % 
% na=numel(qi(1).Attributes);
% for i=1:na
%     qj=qi(1).Attributes(i);
% %     qj.Value
%     if isa(qj.Value,'numeric')
%         valstr=num2str(qj.Value');
%     else
%         valstr=qj.Value.Data;
%     end;
%     disp([qj.Shortname '  ' valstr]);
% end;

