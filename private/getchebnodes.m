function x = getchebnodes(allnx)
% Generate Chebyshev nodes for barypdstepcb
%    (internal function)

d = uint8(length(allnx));
x = zeros(sum(allnx),1);

aid = uint32(0);
for k = 1:d
  x(1+aid:aid+allnx(k)) = 0.5 - cos(linspace(0,1,allnx(k)) * pi) / 2;
	aid = aid + allnx(k);
end
