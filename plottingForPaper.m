% From https://www.youtube.com/watch?v = Lw4IaddGDRg
figure;
surf(peaks);
newmap = cbrewer('seq','Blues',8);
colormap(newmap)

%%
brighten(0.2); % Try to run section multiple times

%%
figure;
colormap(redgreencmap);
cdata = rand(10,100);
h = imagesc(cdata);
colorbar
