 function yaxis_pi(varargin)
%function yaxis_pi(varargin)
% label y axis with various forms of "pi"
% the argument can be a string with p's in it, or fractions of pi:
% [0 1/2 1] or '0 p/2 p' -> [0 pi/2 pi]
% [-1 0 1] or '-p 0 p' -> [-pi 0 pi]
% etc.
% Jeff Fessler

if length(varargin) == 0
	ticks = '0 p';
elseif length(varargin) == 1
	ticks = varargin{1};
else
	error 'only one arg allowed'
end

if ischar(ticks)
	str = ticks;
	str = strrep(str, ' ', ' | ');
	str = strrep(str, '*', '');	% we don't need the "*" in label
	ticks = strrep(ticks, '2p', '2*p');
	ticks = strrep(ticks, '3p', '3*p');
	ticks = strrep(ticks, '5p', '5*p');
	ticks = strrep(ticks, 'p', 'pi');
	ticks = eval(['[' ticks ']']);

else

	if same(ticks, [0])
		str = '0';
	elseif same(ticks, [0 1])
		str = '0 | p';
	elseif same(ticks, [0 1/2 1])
		str = '0 | p/2 | p';
	elseif same(ticks, [-1 0 1])
		str = '-p | 0 | p';
	elseif same(ticks, [0 1 2])
		str = '0 | p | 2p';
	else
		error 'this ticks not done'
	end

end

% here is the main part
axisy(min(ticks), max(ticks))
ytick(ticks)
set(gca, 'yticklabel', str, 'fontname', 'symbol')

function is = same(x,y)
if length(x) ~= length(y)
	is = 0;
	return
end
is = all(x == y);
