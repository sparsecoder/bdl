function plotDicts(Dcell, filename)
N = length(Dcell);
W = ceil(sqrt(N));
REM = W^2 - N > W;
H = REM*(W-1)  +  ~REM*W;

fig = figure('visible', 'off');
for i=1:N
    n = size(Dcell{i},2);
    w = ceil(sqrt(n));
    rem = w^2 - n;
    h = (rem > w)*(w-1)  +  (rem <= w)*w;
    
    p = size(Dcell{i},1);
    
    subaxis(W,H,i, 'SpacingVert', 0.04,'SpacingHoriz', 0, 'Padding', 0, 'MarginTop', 0.027, 'MarginBottom', 0, 'MarginRight', 0, 'MarginLeft', 0);
    d = displayPatches([Dcell{i} zeros(p,rem)]);
    imshow(d(1:h*(sqrt(p)+1) + 1, :), 'Border', 'tight');
    title(i, 'FontSize', 6); set(gca,'xtick',[],'ytick',[]);
end
saveas(fig, filename, 'pdf')