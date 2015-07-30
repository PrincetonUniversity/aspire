function [validx,str]=two_options_question(prompt,val1,val2,defval,valfmt)

val=[];
prompt=sprintf(['%s(' valfmt '/' valfmt ')'],prompt,val1,val2);
if ~isempty(defval)
    prompt=sprintf(['%s[' valfmt '] '],prompt,defval);
end

val1str=sprintf(valfmt,val1);
val2str=sprintf(valfmt,val2);

while isempty(val)
    str=input(prompt,'s');
    if isempty(str)
        val=defval;
    else
        val=sscanf(str,valfmt);
    end
    
    
    valstr=sprintf(valfmt,val);
    
    if strcmp(valstr,val1str)==1
        validx=1;
    elseif strcmp(valstr,val2str)==1
        validx=2;
    else
        val=[];
    end
end