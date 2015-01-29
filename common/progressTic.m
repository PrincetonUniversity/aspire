function progressTic(j,N)
if mod(j,round(N/100))==0
    fprintf('\b|\n');
end