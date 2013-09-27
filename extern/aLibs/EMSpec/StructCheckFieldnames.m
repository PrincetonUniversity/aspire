function s=StructCheckFieldnames(names,s)

if nargin<1
    s=struct;
end;

% check for unknown fields
fnames=fieldnames(s);
for i=1:numel(fnames)
    fname=fnames{i};
    if ~any(strcmp(fname,names))
        warning(['Unknown option ' fname]);
    end;
end;
