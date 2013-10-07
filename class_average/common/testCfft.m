% function testCfft
%
% Tests the functions cfft, cfrft, icfft.
% See also cfft, cfrft, icfft.
%
% Yoel Shkolnisky 9/2/02

function testCfft

eps = 0.00000001;
for i=[100 101 256 317],
	disp(strcat('testing with i=',num2str(i)));   
	x = rand(1,i);

	cdftx  = cdft(x);
	cfftx  = cfft(x);
	frftx  = frft(x,1);
	cfrftx = cfrft(x,1);

	%check frft all this vectors should be equal.

	if isempty(find(cdftx-cfftx>eps))
	     disp('cfft OK');
	else disp('cfft NOT OK');
	end

	if isempty(find(cdftx-frftx>eps))
	     disp('frft OK');
	else disp('frft NOT OK');
	end

	if isempty(find(frftx-cfrftx>eps))
	     disp('cfrft OK');
	else disp('cfrft NOT OK');
	end

	icfftx = icfft(x);
	icdftx = icdft(x);

	if isempty(find(icdftx-icfftx>eps))
	     disp('icfft OK');
	else disp('icfft NOT OK');
   end
   
   disp('**********************');
end
 

    

