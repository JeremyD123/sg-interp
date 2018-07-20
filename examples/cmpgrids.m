function cmpgrids(n,d)
% CMPGRIDS   Compare the available sparse grid types
%    CMPGRIDS(N, D) Compares the maximum-norm-based grid, the
%    no-boundary-nodes grid, the Clenshaw-Curtis grid, the Chebyshev
%    grid and the Gauss-Patterson grid in d = 2 or d = 3 dimensions,
%    level N. If N is omitted, N is taken to 4. If D is omitted, D
%    is taken to 2.
	
if nargin < 1, n = 4; end
if nargin < 2, d = 2; end

set(gcf,'Color',[1 1 1]);
subplot(2,3,1);
options = spset('GridType','Maximum');
plotgrid(n,d,options);
title('Maximum-norm-based grid');
xlabel(['Points: ' num2str(spdim(n,d,options)), ', Level: ' ...
				num2str(n)], 'HorizontalAlignment', 'center'); 
	
subplot(2,3,2);
options = spset('GridType','NoBoundary');
plotgrid(n,d,options);
title('No-boundary-nodes grid');
xlabel(['Points: ' num2str(spdim(n,d,options)), ', Level: ' ...
				num2str(n)], 'HorizontalAlignment', 'center'); 
subplot(2,3,3);
plotgrid(n,d);
title('Clenshaw-Curtis grid');
xlabel(['Points: ' num2str(spdim(n,d)), ', Level: ' ...
				num2str(n)], 'HorizontalAlignment', 'center'); 
subplot(2,3,4);
options = spset('GridType','Chebyshev');
plotgrid(n,d,options);
title('Chebyshev-Gauss-Lobatto grid');
xlabel(['Points: ' num2str(spdim(n,d)), ', Level: ' ...
				num2str(n)], 'HorizontalAlignment', 'center'); 
subplot(2,3,5);
options = spset('GridType','Gauss-Patterson');
plotgrid(n,d,options);
title('Gauss-Patterson grid');
xlabel(['Points: ' num2str(spdim(n,d,options)), ', Level: ' ...
				num2str(n)], 'HorizontalAlignment', 'center'); 
