function cmap = GrabColormap(filenm)

  if exist(filenm)
    cmap = load(filenm);
    cmap = cmap.hot_desat;
  else
    disp(sprintf('Colormap does not exist. Available custom colormaps are: \n'));
    dirout = dir([getenv('FLASHAVM') '/ColorMaps']);
    for i = 1:length(dirout)
      if contains(dirout(i).name, '.mat')
        disp(sprintf('%s\n', dirout(i).name));
      end
    end
    error('');
  end

end
