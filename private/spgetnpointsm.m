function [totalpoints, npoints] = spgetnpointsm(levelseq);
% SPGETNPOINTSM   Compute the number of subgrid points (Max-Grid)
%    (internal function)

nlevels = uint32(size(levelseq, 1));
d = uint16(size(levelseq, 2));
npoints = zeros(nlevels,1,'uint32');

totalpoints = uint32(0);
% Compute the number of points per subdomain of new levels
for kl = 1:nlevels
	np = uint32(1);
	lval = uint8(0);
	for k = 1:d
		lval = levelseq(kl,k);
		if lval == 0 
			np = np * 3;
		else
			np = np * 2^uint32(lval);
		end
	end
	npoints(kl) = np;
	totalpoints = totalpoints + np;
end	
