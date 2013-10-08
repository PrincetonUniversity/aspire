% function testChirp
%
% Tests the functions slowChirp and ChirpZ
%
% Yoel Shkolnisky 13/1/03

function testChirp

% Compare slowChirp to odd fft.
% Compute the aliased Fourier transform for an odd sequence using slowChirp and verify the result.
fprintf('Test 1: SlowChirp odd \t\t');
x = rand(1,7);
A = 1;
W = exp(-2*pi*i/length(x));
c1 = cfft(x);
c2 = slowChirp(x,A,W,length(x));
fprintf('%s\n',assert(c1,c2));

% Compare slowChirp to even fft
% Compute the aliased Fourier transform for an even sequence using slowChirp and verify the result.
fprintf('Test 2: SlowChirp even \t\t');
x = rand(1,8);
A = 1;
W = exp(-2*pi*i/length(x));
c1 = cfft(x);
c2 = slowChirp(x,A,W,length(x));
fprintf('%s\n',assert(c1,c2));

% After this point we assume that slowChirp work correctly.
% In the following tests we use slowChirp as a reference.

% Compare slowChirp and ChirpZ for odd random vector.
% Input length is equal to output length
fprintf('Test 3: ChirpZ odd \t\t');
x = rand(1,11);
A = exp(i*6.5);
W = exp(-2.31*pi*i);
c1 = slowChirp(x,A,W);
c2 = ChirpZ(x,A,W);
fprintf('%s\n',assert(c1,c2));

% Compare slowChirp and ChirpZ for even radon vector
% Input length is equal to output length
fprintf('Test 3: ChirpZ even \t\t');
x = rand(1,16);
A = 2;
W = exp(-2.31*pi);
c1 = slowChirp(x,A,W);
c2 = ChirpZ(x,A,W);
fprintf('%s\n',assert(c1,c2));

% Compare slowChirp and ChirpZ for odd random vector.
% Input length is different from output length
fprintf('Test 4: ChirpZ odd (n<>M) \t');
x = rand(1,21);
A = exp(i*0.1);
W = exp(-2.31*i);
c1 = slowChirp(x,A,W,51);
c2 = ChirpZ(x,A,W,51);
fprintf('%s\n',assert(c1,c2));

% Compare slowChirp and ChirpZ for even random vector.
% Input length is different from output length
fprintf('Test 5: ChirpZ even (n<>M) \t');
x = rand(1,32);
A = exp(i*0.1);
W = exp(-2.31*i);
c1 = slowChirp(x,A,W,64);
c2 = ChirpZ(x,A,W,64);
fprintf('%s\n',assert(c1,c2));


function res = assert(c1,c2)
diff=abs(c1-c2);
if max(diff(:)) > 1.e-12
   res = 'FAILED';
else
   res = 'OK';
end
   
