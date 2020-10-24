clc, clear all

% hawley zabusky test
%db1 = '/data1/blg13/SimulationData/HawleyZabusky/data/n128y/ppm-n128y-uni/hz_hdf5_chk_';
%db2 = '/data1/blg13/SimulationData/HawleyZabusky/data/n128y/mrppm-n128y-uni-kap0e00/hz_hdf5_chk_';
db1 = '/data2/blg13/SimulationData/WDStir/data/n1024/cori/Z01024ppmc_eps0e00_2d_C0.5e0.70e15_fixedEint_tburn_d1e7te1e9edt1e-1cfl0.6/tburn_hdf5_plt_cnt_';
db2 = '/data2/blg13/SimulationData/WDStir/data/n1024/cori/Z01024ppmc_eps1e04_2d_C0.5e0.70e15_fixedEint_tburn_d1e7te1e9edt1e-1cfl0.6/tburn_hdf5_plt_cnt_';

% blast2 test
%db1 = '/data1/blg13/SimulationData/Blast2/data/n512/ppm-n512-uni/blast2_hdf5_chk_';
%db2 = '/data1/blg13/SimulationData/Blast2/data/n512/mrppm-n512-uni-kap1e02/blast2_hdf5_chk_';

% set variable of interest
var = 'vort';

% loop through all checkpoint files
cnt = 0; search = true;
while(search)

  % current files 1 and 2
  cfile1 = sprintf('%s%0.4i',db1,cnt);
  cfile2 = sprintf('%s%0.4i',db2,cnt);

  % compute error between databases (get handles)
  H1 = hdf5info(cfile1);
  H2 = hdf5info(cfile2);

  % get metadata
  meta = GrabHDF5(cfile1, {'node type'});
  nodetyp = meta{1};

  % find desired dataset index
  for i = 1:length(H1.GroupHierarchy.Datasets)
    if strcmp(H1.GroupHierarchy.Datasets(i).Name, sprintf('/%s',var))
      ind = i;
    end
  end

  % read data
  data1 = hdf5read(H1.GroupHierarchy.Datasets(ind));

  % find desired dataset index
  for i = 1:length(H2.GroupHierarchy.Datasets)
    if strcmp(H2.GroupHierarchy.Datasets(i).Name, sprintf('/%s',var))
      ind = i;
    end
  end

  % read data
  data2 = hdf5read(H2.GroupHierarchy.Datasets(ind));

  % check if array sizes are equal
  if (size(data1) == size(data2))

    % compute difference
    %data3 = abs(data1-data2) ./ abs(data1);
    data3 = abs(data1-data2);

    % leaf block data
    leaf = find(nodetyp==1);
    leafdata = data3(:,:,1,leaf);

    % compute norms
    L1(cnt+1) = norm(leafdata(:), 1);
    L2(cnt+1) = norm(leafdata(:), 2) / sqrt(length(leafdata(:)));
    Linf(cnt+1) = norm(leafdata(:), 'inf');
    Linf(cnt+1) = max(abs(leafdata(:)));

  else

    % error message
    disp('Meshes are not equal.');

    % just set to zeros
    data3 = zeros(size(data1));

  end

  % append resultant data to one of the hdf5 files
  h5write(cfile1,'/velz',data3);
  %h5write(cfile1,'/shck',data3);

  % check if file dne
  if exist(sprintf('%s%0.4i',db1,cnt+1)) ...
    + exist(sprintf('%s%0.4i',db2,cnt+1)) == 0

    % stopping
    search = false;

  end

  % update counter
  cnt = cnt + 1;

end

% plot the L2 error
figure(1)
%semilogy(Linf,'k');
%semilogy(L2,'k'); hold on;

%x = 
