if ~is_gpu()
    warning('No GPU is present. Skipping compilation of CUDA MEX files.');
    return;
end

debug=1;
timing=0;

flags='-O -g -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64/ -lcublas -lcudart';

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
