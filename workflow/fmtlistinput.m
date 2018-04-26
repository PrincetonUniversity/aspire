function [vals,str]=fmtlistinput(prompt,defval,outfmt)
% Input comma separated list
if ~isempty(defval)
    prompt=sprintf(['%s[' outfmt ']'],prompt,defval);
end

vals={};
while isempty(vals) || isempty(vals{1})
    str=input(prompt,'s');
    C=strsplit(str,',');
    vals{1}=defval;
    
    if ~isempty(str)
        if ~isempty(C{1})
            vals=cell(numel(C),1);
            for k=1:numel(C)
                vals{k}=sscanf(C{k},outfmt);
            end
        else % Input was malformed. Ask for input again
            vals=[];
        end
    else % Input was empty. Return default value
        str=sprintf(outfmt,defval);        
    end
end
