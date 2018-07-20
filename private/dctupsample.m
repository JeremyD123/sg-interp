function ip = dctupsample(z, nnew)
% DCTUPSAMPLE  Perform upsampling of 1-dimensional grid data (at 
%    Chebyshev-Gauss-Lobatto nodes) using FFT/IFFT. The function
%    returns only the even nodal values, as required by the
%    hierarchical sparse grid construction algorithm.
%    (internal function)

n = size(z,1);
m = size(z,2);

% Extend grid values to grid points of trigonometric polynomial
zext = [z; z(n-1:-1:2,:)];

% Perform FFT
chebc = fft(zext);

% Rescale weights & pad with zeros
chebc = (chebc(1:n,:)) / (n-1) * (nnew - 1);
chebc(n,:) = chebc(n,:) / 2;
paddchebc = [chebc; zeros(2*(nnew-n)-1,m); chebc(n:-1:2,:)];

% Perform inverse FFT
ipext = real(ifft(paddchebc));

% Remove irrelevant nodes; these are always all the odd indices;
% keep the even ones.
ip = ipext(2:2:nnew,:);
