function [maxMach simtime] = ComputeMaxMachNo(filenm)

  vars = {'velx','dens','pres','gamc'};

  % construct cell array to pass
  inputs = [{'node type', 'bounding box', 'refine level', 'real scalars'} vars];

  % get 'varnm' variable data from current hdf5 data
  tmpData = GrabHDF5(filenm, inputs);

  % number of meta data elements
  nmeta = 4;

  % get simulation time
  for m = 1:length(tmpData{4})
    scalarnm = strtrim(h5stringconvert(tmpData{4}(m).Data{1}));
    if strcmp(scalarnm, 'time')
      simtime = tmpData{4}(m).Data{2};
    elseif strcmp(scalarnm, 'dt')
      simdt = tmpData{4}(m).Data{2};
    end
  end


  % get some mesh metadata
  nodetyp = tmpData{1};                  % the type of AMR block (leaf or not)
  nblocks = length(tmpData{1});          % the number of AMR blocks in the mesh
  lrefine = tmpData{3};                  % the refinement level of the block ID

  % loop through all blocks
  maxMach = 0.0;
  for blk = 1:nblocks

    % check if leaf block
    if nodetyp(blk) == 1

      velx = tmpData{nmeta+1}(:,1,1,blk);
      dens = tmpData{nmeta+2}(:,1,1,blk);
      pres = tmpData{nmeta+3}(:,1,1,blk);
      gamc = tmpData{nmeta+4}(:,1,1,blk);

      c0 = sqrt(gamc.*pres./dens);
      maxMach = max(maxMach, max(abs(velx) ./ c0));

    end

  end

end
