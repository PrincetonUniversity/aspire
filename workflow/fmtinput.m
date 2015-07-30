function [val,str]=fmtinput(prompt,defval,outfmt)
if ~isempty(defval)
    prompt=sprintf(['%s[' outfmt ']'],prompt,defval);
end

val=[];
while isempty(val)
    str=input(prompt,'s');
    val=defval;

    if ~isempty(str)
        val=sscanf(str,outfmt);
    else % Input was empty. Return default value
        str=sprintf(outfmt,defval);
    end
end
