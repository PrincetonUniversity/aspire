debug=0;
timing=0;

flags='-O -I/usr/local/cuda/targets/x86_64-linux/include/ -L/usr/local/cuda/targets/x86_64-linux/lib/ -lcublas';

debugdef='';
if debug
debugdef= '-DDEBUG';
end

timingdef='';
if timing
    timingdef='-DTIMING';
end

files=dir('gpuaux*.cpp');

for k=1:numel(files)
    fname=files(k).name;
    cmd =['mex ' fname ' ' flags ' ' debugdef ' ' timingdef];
    eval(cmd);
end
