function [ipder] = ppderiv(ipder, ipder2, maxlev, y)
% PPDERIV  Derivative post processing. Combine two computed 
%    derivatives to one via linear interpolation.
%    (internal function)

d = size(y,2);
ninterp = uint32(size(y,1));

if length(maxlev) < d
  % Expand maxLevel to a vector in case we have a
	% conventional, non-adaptive sparse grid.
	maxlev = ones(d,1,'uint8') * maxlev;
end

for k = 1:ninterp
	for l = 1:d
	  if maxlev(l) == 0, continue; end
		yt = y(k,l);
		stepsize = 1/2^double(maxlev(l));
		halfstep = 0.5 * stepsize;
		if maxlev(l) == 1
			ytd1 = 0.25;
		elseif maxlev(l) > 1
			if yt <= halfstep
				ytd1 = halfstep;
			elseif yt >= 1 - halfstep;
				ytd1 = 1 - 1.5*stepsize;
			else																							                     
				ytd1 = halfstep + ...
				       floor( (yt - halfstep) / stepsize) * stepsize;
			end
		end
		ipd1 = ipder(k,l);
		ipd2 = ipder2(k,l);
		if ipd1 * ipd2 >= 0 
			ipder(k,l) = ipd1 + (ipd2 - ipd1) / stepsize * (yt - ytd1);
		elseif yt <= ytd1 + halfstep
			ipder(k,l) = ipd1 - ipd1 / halfstep * (yt - ytd1);
		else
			ipder(k,l) = ipd2 / halfstep * (yt - ytd1 - halfstep);
		end
	end
end
