function [totalpoints, npoints] = spgetnpointsnb(levelseq);
% SPGETNPOINTSNB   Compute the number of subgrid points (NB-Grid)
%    (internal function)

nlevels = uint32(size(levelseq, 1));
d = uint16(size(levelseq, 2));
npoints = zeros(nlevels,1,'uint32');

totalpoints = uint32(0);
% Compute the number of points per subdomain of new levels
for kl = 1:nlevels
	np = uint32(1);
	for k = 1:d
		np = np * 2^uint32(levelseq(kl,k));
	end
	npoints(kl) = np;
	totalpoints = totalpoints + np;
end	
