path = 'stem/5nm/';
psize = 24^2;
p.n1 = psize^0.5;
p.n2 = psize^0.5;
p.d1 = 2;
p.d2 = 2;

dl = {};
for i=1:18
    fname = sprintf('%s%d.tif',path,i);
    fprintf('%s',fname);
    img = Image(fname, p);
    dl{i} = BPFA(img.patches, img.patches, 9);
    clear img
    dl{i}.learn(20);
end

