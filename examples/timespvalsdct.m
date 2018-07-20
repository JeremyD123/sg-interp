function timespvalsdct
% TIMESPVALSDCT   Compares the performance of SPVALS for the Cheby-
%    shev-Gauss-Lobatto grid with and without the DCT-based grid
%    construction algorithm.
%
% Note: This demo takes a couple of minutes to run.
%
%    See also TIMESPVALS.

n1 = [14,12,10,6];
n2 = [17,14,10,6];
dvec = [1,2,4,8];

for k = 1:length(n1)
  subplot(2,2,k,'align');
  d = dvec(k);
  options = spset('GridType', 'Chebyshev', 'SparseIndices', 'on', ...
								'EnableDCT', 'off');
  timespvals(options,n1(k),d);
  options = spset('GridType', 'Chebyshev', 'SparseIndices', 'on');
  hold on
  timespvals(options,n2(k),d);
  hold on
  plot([1,1e7],[1e-4,1e3],'k-.')
  hold off
  title(['d = ' num2str(d)]);
  xlabel('');
  axis normal;
  h = get(gca,'Children');
  set(h(2),'LineStyle','--','LineWidth',1.5,'Marker','*','Color','r');
  set(h(3),'LineStyle','-','LineWidth',1.5,'Marker','+','Color','b');
  legend({'no DCT', 'with DCT', 'O(N)'},2);
  set(gca,'XTick',[1e1,1e2,1e3,1e4,1e5]);
  set(gca,'YTick',[1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4]);
  set(gca,'XLim',[1e1,1e5])
end
hold off;
