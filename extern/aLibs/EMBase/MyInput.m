function val=MyInput(text,val)
% Input a numeric value.
% Returns the original value if nothing was typed.
% This function can be used to accept a row vector too.

str=input([text ' [' num2str(val) ']? '],'s');
if numel(str)==0
    return
else
    newVal=str2num(str);
    if numel(newVal)>0
        val=newVal;
    end;
end;
