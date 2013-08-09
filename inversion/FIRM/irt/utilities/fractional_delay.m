 function y = fractional_delay(x, delay)
%function y = fractional_delay(x, delay)
%
% given N samples x[n] of a real, periodic, band-limited signal x(t),
% compute sinc interpolated samples of delayed signal y(t) = x(t - delay)
% each column of x can be shifted by a different amount if delay is a vector.
% in
%	x	[N,L]
%	delay	[L]
% out
%	y	[N,L]
%
% see laakso:96:stu for more ideas.
%
% Copyright 2003-11-1, Jeff Fessler, The University of Michigan
% Extend to allow x to have multiple columns, 2003-11-2, Yingying Zhang.

if nargin < 1, help(mfilename), disp('Need Inputs'), error(mfilename), end
if streq(x, 'test'), fractional_delay_test, return, end
if size(x,2) ~= length(delay)
	error 'Need size(x) = [N,L]; length(delay) = L'
end

dims = size(x);
N = dims(1);
L = dims(2);
if length(dims) > 2, error 'x must be 1d or 2d', end
X = fft(x);	% fft of each column

% it is important to choose the k indices appropriately!
if rem(N,2)	% odd
	k = [-(N-1)/2:(N-1)/2]';
else
	k = [-N/2:(N/2-1)]';
end

c = exp(-1i * 2*pi/N * k * delay(:)');	% [N,L] outer product

if ~rem(N,2)	% even
	mid = 1;
	c(mid,:) = real(c(mid,:));	% this is the other key trick!
end

c = ifftshift1(c);	% ifftshift differs from fftshift for odd N!

Y = X .* c;
y = ifft(Y);


%
% 1D ifftshift for each column
%
function c = ifftshift1(c)
[N,L] = size(c);
if L == 1
	c = ifftshift(c);	% ifftshift differs from fftshift for odd N!
else % multiple input signals
	if ~rem(N,2)	% even
		c = c([N/2+1:end 1:N/2],:);
	else		% odd
		c = c([(N-1)/2+1:end 1:(N-1)/2],:);
	end
end


%
% self test
%
function fractional_delay_test

Nlist = [5 6];
im clf, pl=240;
for ii=1:2
	N = Nlist(ii);
	n = [0:(N-1)]';
	xt = inline('sinc_periodic(t, N)', 't', 'N');
	x = xt(n, N);
	xx = [x x];
	delay = [3.7; -2.2];
	y = fractional_delay(xx, delay);

	t = linspace(0,2*N,401);
	yt1 = xt(t-delay(1),N);
	yt2 = xt(t-delay(2),N);

	if im
	subplot(pl+0+ii), plot(t, xt(t,N), '-', n, xx(:,1), 'o')
	axis([0 2*N -0.4 1.1]), title(sprintf('N=%d 1st input', N))
	subplot(pl+2+ii), plot(t, xt(t,N), '-', n, xx(:,2), 'o')
	axis([0 2*N -0.4 1.1]), title(sprintf('N=%d 2nd input', N))
	subplot(pl+4+ii)
	plot(t, yt1, '-', n, real(y(:,1)), 's', n, imag(y(:,1)), '.')
	axis([0 2*N -0.4 1.1]), title(sprintf('delay=%g 1st input',delay(1)))
	subplot(pl+6+ii)
	plot(t, yt2, '-', n, real(y(:,2)), 's', n, imag(y(:,2)), '.')
	axis([0 2*N -0.4 1.1]), title(sprintf('delay=%g 2nd input',delay(2)))
	end
end
