function plotDicts(Dcell, filename)
N = length(Dcell);
W = ceil(sqrt(N));
H = height(W,N);

fig = figure('visible', 'off');
for i=1:N
    n = size(Dcell{i},2);
    w = ceil(sqrt(n));
    [h, dif] = height(w,n);
    
    p = size(Dcell{i},1);
    
    subaxis(W,H,i, 'SpacingVert', 0.04,'SpacingHoriz', 0, 'Padding', 0, 'MarginTop', 0.027, 'MarginBottom', 0, 'MarginRight', 0, 'MarginLeft', 0);
    d = displayPatches([Dcell{i} zeros(p,dif)]);
    imshow(d(1:h*(sqrt(p)+1) + 1, :), 'Border', 'tight');
    title(i, 'FontSize', 6);
    set(gca,'xtick',[],'ytick',[]);
end

saveas(fig, filename, 'pdf')

function [h, dif] = height(w,n)
dif = w^2 - n;
h = (dif > w)*(w-1)  +  (dif <= w)*w;
