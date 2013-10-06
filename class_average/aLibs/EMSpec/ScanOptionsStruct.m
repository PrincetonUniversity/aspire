function options=ScanOptionsStruct(options,validOptions,optionDefaults)
% function opts=ScanOptionsStruct(options,validOptions,optionDefaults)
% Scan the fields of the struct options.  If there is a field that doesn't
% match a member of the validOptions cell array, issue a warning.  If there
% validOptions{i} is a name that is not represented as an options field,
% create the field and set its value to optionDefaults(i), or zero if the
% third argument isn't given.

nv=numel(validOptions);
if nargin<3
    optionDefaults=zeros(nv,1);
end;

optNames=fieldnames(options);
nopt=numel(optNames);
for i=1:nopt
    matches=strcmp(optNames{i},validOptions);
    if ~any(matches)
        warning(['Invalid option: ' optNames{i}]);
    end;
end;
for i=1:nv
    matches=strcmp(validOptions{i},optNames);
    if ~any(matches)
        options=setfield(options,validOptions{i},optionDefaults(i));
    end;
end;
