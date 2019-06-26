
function [suff]=datetime_suffix()
    s=split(datestr(datetime),' ');
    suff=strcat(strrep(s(1),'-',''),'_',strrep(s(2),':',''));
end