 function ob = power(ob, sup)
%function ob = power(ob, sup)
% array "power" method (G.^2)
% matrix power (mpower) would be for G^2

if isempty(ob.handle_power)
	error 'no power method for this object'
end

ob = feval(ob.handle_power, ob, sup);
if isfield(ob, 'scale') && ob.scale ~= 1
	ob.scale = ob.scale .^ sup;
end
