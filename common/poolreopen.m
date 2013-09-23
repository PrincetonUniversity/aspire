function poolreopen(n)

if nargin<1
    n=8;
end

ps=matlabpool('size');
if ps==0 || ps~=n
    if ps~=0
        matlabpool close;
    else
        matlabpool('open',n)
    end
end
