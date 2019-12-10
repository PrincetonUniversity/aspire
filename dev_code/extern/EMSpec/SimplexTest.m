% Example of use:  minimize (p-q).^6
p=[2 2 1.5 2]';  % inital guess
q=[1 1 1 1]';       % true values
[p,t]=Simplex('init',p);
iters=0;
for i=1:200
    y=sum((p-q).^6);
    [p,t]=Simplex(y);
    if ~mod(i, 10)  % every 10 iterations print out the fitted value
        pc=Simplex('centroid')'
    end;
end;
p=Simplex('centroid');
