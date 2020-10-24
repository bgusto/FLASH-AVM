function Data = GrabHDF5(filename, vars)

  % read data
  hinfo = hdf5info(filename);

  % find desired datasets
  for i = 1 : length(hinfo.GroupHierarchy.Datasets)
    for j = 1 : length(vars)
      if strcmp( hinfo.GroupHierarchy.Datasets(i).Name, sprintf('/%s',vars{j}) )
        Data{j} = hdf5read( hinfo.GroupHierarchy.Datasets(i) );
      end
    end
  end

end
