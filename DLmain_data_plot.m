path = '/Users/donguk.kim/projects/coupling_xolotl/';

colIDforx = 1;
colIDfory = 2;

nprocs = 4;

for i=0:nprocs-1
    fname = sprintf('dataout_%02d.dat',i);
    path_fname = [path fname];
    dataval = dlmread(path_fname);
    nrow = size(dataval,1);
    ncol = size(dataval,2);
    xdat = dataval(1:nrow, colIDforx);
    ydat = dataval(1:nrow, colIDfory);
    if (i == 0)
        xall = xdat;
        yall = ydat;
    else 
        xall = cat(1, xall, xdat);
        yall = cat(1, yall, ydat);
    end
end

xmin = min(xall);
xmax = max(xall);
ymin = min(yall);
ymax = max(yall);

plot(xall,yall, '-k', 'linewidth',1.5);
xlabel('x', 'fontsize', 22);
ylabel('y', 'fontsize', 22);
set(gca,'fontsize', 22);
xlim([xmin xmax]);
ylim([ymin ymax]);
