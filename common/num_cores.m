%CORES Count the number of cores available

%

% Usage

%    n = num_cores();

%
% Output
%    n: Number of cores available for parallel computing.

function n = num_cores()
	v = ver('distcomp');
	if isempty(v.Version)
		n = 1;
	elseif str2double(v.Version) >= 6.3
		p = gcp;
		n = p.NumWorkers;
	else
		n = matlabpool('size');
	end
	n = max(1, n);
end


