function fullseq = spseq2full(seq)
% SPSEQ2FULL   Convert sparse sequence of index sets to full one
% (internal function)

d = uint16(size(seq.forwardNeighbors, 2));
nlevels = uint32(size(seq.indicesNDims, 1));

fullseq = zeros(nlevels, d, 'uint8');
addr = uint32(1);
for k = 1:nlevels
	ndims = seq.indicesNDims(k);
	for did = 1:ndims
		fullseq(k,seq.indicesDims(addr)) = seq.indicesLevs(addr);
		addr = addr + 1;
	end
end
