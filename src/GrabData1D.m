function [vardata x xc delx] = GrabData1D(filenm, varnm, intrp, reflvl)
%
% GrabData1D Grab one-dimensional data from the specified input file.
%
%-------------------------------------------------------------------------------%
% Info:
%   This function grabs data from FLASH hdf5 files for the variable of interest.
%   The function takes block-partitioned data and converts it into a single
%   array. For multi-level adaptive mesh refinement data, the user may choose to
%   prolong all of the data to the finest level using the 'intrp' argument. The
%   function assumes that the blocks contain a fixed number of cells. Currenty,
%   the only choice for prolonging data is a direct copy. Third-order
%   average-preserving interpolation coming soon.
%
% Inputs:
%   filenm - the hdf5 filename
%   varnm  - the variable name to be plotted from hdf5 files (a string)
%   intrp  - optional interpolation choice, if prolonging data
%   reflvl - the desired refinement level, using the coarsest parent block
%            found in the input file as the base level (found automatically if
%            not provided and compatible 'intrp' option specified)
%
% Outputs:
%   vardata - the output array
%   x       - the cell edges
%   xc      - the cell centers
%   dx      - the cell widths
%
%-------------------------------------------------------------------------------%

  % determine if we want to prolong AMR data to finest level
  if nargin < 3
    prolong = false;
  else
    prolong = true;
  end

  % get 'varnm' variable data from current hdf5 data
  Data = GrabHDF5(filenm,{'node type', 'bounding box', 'refine level', varnm});

  % mesh metadata 
  nodetyp = Data{1};                  % the type of AMR block (leaf or not)
  nblocks = length(Data{1});          % the number of AMR blocks in the mesh
  lrefine = Data{3};                  % the refinement level of the block ID
  nxb = length(Data{4}(:,1,1,1));     % the number of cells per block

  % initialize block bounding coordinates array
  mdata = zeros(2,length(find(nodetyp==1)));

  % determine max level automatically if prolonging and no desired refinement level given
  if prolong
    if nargin < 4
      reflvl = max(lrefine);
    end
  end

  % initialize output arrays for the non-prolonged case
  if ~prolong
    x = [];
    xc = [];
    delx = [];
    vardata = [];
  end

  % initialize the global xmin, xmax
  xmin = 0.0; xmax = 0.0;
  maxlvl = 0;

  % initialize minimum dx (because we want to know what dx is on finest level)
  dxmin = 1e25;

  % loop through blocks
  cnt = 1;
  for blk = 1:nblocks

    % check if leaf block
    if nodetyp(blk) == 1

      % ends of block (1=lo,2=hi)
      xlo = Data{2}(1,1,blk);
      xhi = Data{2}(2,1,blk);

      % determine the bounding box of the global block
      xmin = min(xmin,xlo);
      xmax = max(xmax,xhi);

      % cell width
      dx = (xhi - xlo) / nxb;
      dxmin = min(dxmin, dx);

      % remember which level dxmin found on
      maxlvl = max(maxlvl, lrefine(blk));

      % save in mesh data
      mdata(1,cnt) = xlo;
      mdata(2,cnt) = xhi;

      % if not prolonging, concatenate current block information to global adaptive array
      if ~prolong
        x = [x xlo:dx:xhi];
        xc = [xc (xlo+dx/2):dx:(xhi-dx/2)];
        delx = [delx dx*ones(1,nxb)];
        vardata = [vardata Data{4}(:,1,1,blk)'];
      end

    end

    % increment counter
    cnt = cnt + 1;

  end

  % generate output data (prolong or not)
  if prolong

    % determine number of base blocks
    nbase = 0;
    for blk = 1:nblocks
      if lrefine(blk) == 1
        nbase = nbase + 1;
      end
    end

    % create the global domain, compute number of cells
    if reflvl == maxlvl

      x = xmin:dxmin:xmax;
      nx = length(x)-1;

    else

      % compute number of cells on finest level to cover domain
      nx = nbase * nxb * 2^(reflvl-1);
      x = linspace(xmin, xmax, nx+1);

    end

    % compute cell centers and cell widths as well
    xc = abs(x(1:end-1) + x(2:end)) / 2;
    delx = abs(x(2:end) - x(1:end-1));

    % initialize the global block of data
    vardata = zeros(1,nx);

    % now map the data on blocks to the global data array
    for blk = 1:nblocks

      % check if leaf block
      if nodetyp(blk) == 1

        % ends of block (1=lo,2=hi)
        xlo = Data{2}(1,1,blk);
        xhi = Data{2}(2,1,blk);

        % get data
        blkdata = Data{4}(:,1,1,blk);

        % get amr level
        lvl = lrefine(blk);

        % prolong the data to finest level
        for l = lvl:reflvl-1

          % get size of current block
          nxblk = length(blkdata);

          % create new data for block at higher level
          blkdata2 = zeros(1,2*nxblk);

          % determine how to prolong data
          if strcmp(intrp,'copy')

            % copy parent cells to children
            for i = 1:nxblk
              blkdata2(2*i-1) = blkdata(i);
              blkdata2(2*i) = blkdata(i);
            end

          else

            error('Invalid prolongation option: intrp. Valid options are ... copy')

          end

          % copy variable for next iteration
          blkdata = blkdata2;

        end

        % get size of this newly created block
        nxblk = length(blkdata);

        % determine the starting (i,j) indices in global array
        iglb = int32((xlo - xmin) / dxmin + 1);

        % copy data
        for i = 0:nxblk-1
          vardata(iglb+i) = blkdata(i+1);
        end

      end

    end

  else

    % sort data based on coordinates
    x = unique(x);
    [xc ia] = sort(xc);
    delx = vardata(ia);
    vardata = vardata(ia);

  end

end
