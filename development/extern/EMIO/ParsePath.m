function [base final]=ParsePath(pathText)
% function [base final]=ParsePath(pathText);
% From a path string, extract the final directory from it.
% For example, ParsePath('/Users/fred/data') returns
% base='/Users/fred/' and final='data/'.  Both returned strings (if not null
% strings) end with a slash.

% defaults
base=AddSlash(pathText);
    final='';

nt=numel(pathText);
p=strfind(pathText,'/');
np=numel(p);
if np==0
    return
end;
if p(np)==nt  % ignore trailing slash
    np=np-1;
end;
if np>0 && p(np)<nt
    final=pathText(p(np)+1:nt);
    final=AddSlash(final);
    base=pathText(1:p(np));
end;
