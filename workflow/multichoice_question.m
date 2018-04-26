function [val,str]=multichoice_question(prompt,options,vals,defstr)

% Convert options to a cell array
if ~iscell(options)
    error('options must be a cell array of strings')
end
    

prompt=sprintf('%s(%s',prompt,options{1});
for k=2:numel(options)
    prompt=sprintf('%s/%s',prompt,options{k});
end
prompt=sprintf('%s)',prompt);
if ~isempty(defstr)
    prompt=sprintf('%s? [%s] ',prompt,defstr);
end

str=[];
while isempty(str)
    str=input(prompt,'s');
    if isempty(str)
        str=defstr;
    end

    matches=strcmp(str,options);
    matchidx=find(matches==1,1,'first');
    
    if isempty(matchidx)
        str=[];
    else
        val=vals(matchidx);
    end
end