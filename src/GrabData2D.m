function [vdata mdata x y] = GrabData2D(filenm, varnm, intrp)
%
% GrabData2D Grab one-dimensional data from the specified input file.
%
%-------------------------------------------------------------------------------%
% Info: This function grabs data from FLASH hdf5 files for the variable of
%   interest.  The function returns a an array the size of the numbre of blocks,
%   with each element containing that block's structured data (specified by
%   'varnm'). Assumes data is cell-centered (no face-vars).
%
% Inputs:
%   filenm -  the hdf5 filename
%   varnm -   the variable name to be plotted from hdf5 files (a string)
%   intrp  - optional interpolation choice, if prolonging data
%
% Outputs:
%   vdata - the output data matrix
%   mdata - bounding box for each leaf block
%   x - meshgrid, x
%   y - meshgrid, y
%
%-------------------------------------------------------------------------------%

  % get 'varnm' variable data from current hdf5 data
  Data = GrabHDF5(filenm,{'node type', 'bounding box', 'refine level', varnm});

  % amr metadata 
  nodetyp = Data{1};
  nblocks = length(Data{1});
  lrefine = Data{3};
  nxb = length(Data{4}(:,1,1,1));
  nyb = length(Data{4}(1,:,1,1));

  % initialize block bounding coordinates array
  mdata = zeros(4,length(find(nodetyp==1)));

  % initialize xmin, xmax, ymin, ymax + dxmin + dymin
  xmin = 1e25; xmax = -1e25;
  ymin = 1e25; ymax = -1e25;

  % initialize to large number
  dxmin = 1e25;
  dymin = 1e25;

  % initial sweep through data to get xmin, xmax, ymin, ymax + dx, dy
  cnt = 1;
  for blk = 1:nblocks

    % check if leaf block
    if nodetyp(blk) == 1

      % ends of block (1=lo,2=hi)
      xlo = Data{2}(1,1,blk);
      xhi = Data{2}(2,1,blk);
      ylo = Data{2}(1,2,blk);
      yhi = Data{2}(2,2,blk); 

      % determine the bounding box of the global block
      xmin = min(xmin,xlo);
      xmax = max(xmax,xhi);
      ymin = min(ymin,ylo);
      ymax = max(ymax,yhi);

      % cell width
      dxmin = min(dxmin, (xhi - xlo) / nxb);
      dymin = min(dymin, (yhi - ylo) / nyb);

      % save in mesh data
      mdata(1,cnt) = xlo;
      mdata(2,cnt) = xhi;
      mdata(3,cnt) = ylo;
      mdata(4,cnt) = yhi;

      % increment counter
      cnt = cnt + 1;

    end

  end

  % create the global domain
  x = xmin:dxmin:xmax;
  y = ymin:dymin:ymax;

  % number of cells
  nx = length(x)-1;
  ny = length(y)-1;

  % create meshgrid
  [x, y] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));

  % initialize the global block
  vdata = zeros(ny,nx);

  % now map the data on blocks to the global data array
  for blk = 1:nblocks

    % check if leaf block
    if nodetyp(blk) == 1

      % ends of block (1=lo,2=hi)
      xlo = Data{2}(1,1,blk);
      xhi = Data{2}(2,1,blk);
      ylo = Data{2}(1,2,blk);
      yhi = Data{2}(2,2,blk); 

      % get data
      blkdata = Data{4}(:,:,1,blk);

      % check amr level
      lvl = lrefine(blk);

      % prolong the data to finest level
      for l = lvl:max(lrefine)-1

        % get size of current block
        [nxblk nyblk] = size(blkdata);

        % create new data for block at higher level
        blkdata2 = zeros(2*nyblk,2*nxblk);

        if strcmp(intrp,'copy')

            % loop through cells in block
            for j = 1:nyblk
              for i = 1:nxblk

                % copy parent cells to children
                blkdata2(2*i-1,2*j-1) = blkdata(i,j);
                blkdata2(2*i-1,2*j) = blkdata(i,j);
                blkdata2(2*i,2*j-1) = blkdata(i,j);
                blkdata2(2*i,2*j) = blkdata(i,j);

              end
            end

        elseif strcmp(intrp,'avginterp')

            % loop through cells that use centered stencil
            for j = 2:nyblk-1
                for i = 2:nxblk-1

                    ioff = 0;
                    joff = 0;

                    u11 = blkdata(i-1+ioff,j-1+joff); u12 = blkdata(i-1+ioff,j+joff); u13 = blkdata(i-1+ioff,j+1+joff);
                    u21 = blkdata(i+ioff,j-1+joff);   u22 = blkdata(i+ioff,j+joff);   u23 = blkdata(i+ioff,j+1+joff);
                    u31 = blkdata(i+1+ioff,j-1+joff); u32 = blkdata(i+1+ioff,j+joff); u33 = blkdata(i+1+ioff,j+1+joff);

                    % predict based on formula
                    qx = -(u32 - u12) / 8;
                    qy = -(u23 - u21) / 8;
                    qxy = (u33 - u31 - u13 + u11) / 64;

                    blkdata2(2*i-1,2*j-1) = u22 + qx + qy + qxy;
                    blkdata2(2*i-1,2*j) = u22 + qx - qy - qxy;
                    blkdata2(2*i,2*j-1) = u22 - qx + qy - qxy;
                    blkdata2(2*i,2*j) = u22 - qx - qy + qxy;

                end
            end

            % loop through cells that use lower-biased
            j = nyblk;
            for i = 2:nxblk-1

                ioff = 0;
                joff = -1;

                u11 = blkdata(i-1+ioff,j-1+joff); u12 = blkdata(i-1+ioff,j+joff); u13 = blkdata(i-1+ioff,j+1+joff);
                u21 = blkdata(i+ioff,j-1+joff);   u22 = blkdata(i+ioff,j+joff);   u23 = blkdata(i+ioff,j+1+joff);
                u31 = blkdata(i+1+ioff,j-1+joff); u32 = blkdata(i+1+ioff,j+joff); u33 = blkdata(i+1+ioff,j+1+joff);

                blkdata2(2*i-1,2*j-1) = (u11 - 4.0*u12 - 5.0*u13 + 8.0*u21 - ...
                                            32.0*u22 - 40.0*u23 - 1.0*u31 + 4.0*u32 + ...
                                            5.0*u33) * (-0.015625);
                blkdata2(2*i-1,2*j) = (u11 - 4.0*u12 + 11.0*u13 + 8.0*u21 ...
                                            - 1.0*u31 + 4.0*u32 - ...
                                            11.0*u33 -32.0*u22 + 88.0*u23) * 0.015625;
                blkdata2(2*i,2*j-1) = (u11 - 4.0*u12 - 5.0*u13 - 8.0*u21 ...
                                            - 1.0*u31 + 4.0*u32 + ...
                                            5.0*u33 + 32.0*u22 + 40.0*u23)  * 0.015625;
                blkdata2(2*i,2*j) =  (u11 - 4.0*u12 + 11.0*u13 - 8.0*u21 ...
                                            - 1.0*u31 + 4.0*u32 - ...
                                            11.0*u33 + 32.0*u22 - 88.0*u23) * (-0.015625);

            end

            % lower edge
            j = 1;
            for i = 2:nxblk-1

                ioff = 0;
                joff = 1;

                u11 = blkdata(i-1+ioff,j-1+joff); u12 = blkdata(i-1+ioff,j+joff); u13 = blkdata(i-1+ioff,j+1+joff);
                u21 = blkdata(i+ioff,j-1+joff);   u22 = blkdata(i+ioff,j+joff);   u23 = blkdata(i+ioff,j+1+joff);
                u31 = blkdata(i+1+ioff,j-1+joff); u32 = blkdata(i+1+ioff,j+joff); u33 = blkdata(i+1+ioff,j+1+joff);

                blkdata2(2*i-1,2*j-1) = (11.0*u11 - 4.0*u12 + u13 + 88.0*u21 - ...
                                32.0*u22 + 8.0*u23 - 11.0*u31 + 4.0*u32 - ...
                                1.0*u33) * (0.015625);

                blkdata2(2*i-1,2*j) = (5.0*u11 + 4.0*u12 - 1.0*u13 + ...
                                40.0*u21 + 32.0*u22 - 8.0*u23 - 5.0*u31 - ...
                                4.0*u32 + u33) * (0.015625);

                blkdata2(2*i,2*j-1) = (11.0*u11 - 4.0*u12 + u13 + ...
                                32.0*u22 - 8.0*u23 - 11.0*u31 + 4.0*u32 - ...
                                1.0*u33 - 88.0*u21)  * (-0.015625);

                blkdata2(2*i,2*j) = (5.0*u11 + 4.0*u12 - 1.0*u13 - ...
                                40.0*u21 - 32.0*u22 + 8.0*u23 - 5.0*u31 - ...
                                4.0*u32 + u33) * (-0.015625);

            end

            % right edge
            i = nxblk;
            for j = 2:nyblk-1

                ioff = -1;
                joff = 0;

                u11 = blkdata(i-1+ioff,j-1+joff); u12 = blkdata(i-1+ioff,j+joff); u13 = blkdata(i-1+ioff,j+1+joff);
                u21 = blkdata(i+ioff,j-1+joff);   u22 = blkdata(i+ioff,j+joff);   u23 = blkdata(i+ioff,j+1+joff);
                u31 = blkdata(i+1+ioff,j-1+joff); u32 = blkdata(i+1+ioff,j+joff); u33 = blkdata(i+1+ioff,j+1+joff);

                blkdata2(2*i-1,2*j-1) = 4*((5*u11)/256 - (5*u12)/64 + (55*u13)/256 + u21/64 - u22/16 + ...
                                (11*u23)/64 - u31/256 + u32/64 - (11*u33)/256);
                blkdata2(2*i-1,2*j) = 4*((11*u11)/256 - (11*u12)/64 + (121*u13)/256 - u21/64 + ...
                                u22/16 - (11*u23)/64 + u31/256 - u32/64 + (11*u33)/256);
                blkdata2(2*i,2*j-1) = 4*((11*u12)/64 - (11*u11)/256 + (55*u13)/256 + u21/64 - u22/16 - ...
                                (5*u23)/64 - u31/256 + u32/64 + (5*u33)/256);
                blkdata2(2*i,2*j) = 4*((5*u12)/64 - (5*u11)/256 + (25*u13)/256 - u21/64 + u22/16 + ...
                                (5*u23)/64 + u31/256 - u32/64 - (5*u33)/256);

            end

            % left edge
            i = 1;
            for j = 2:nyblk-1

                ioff = 1;
                joff = 0;

                u11 = blkdata(i-1+ioff,j-1+joff); u12 = blkdata(i-1+ioff,j+joff); u13 = blkdata(i-1+ioff,j+1+joff);
                u21 = blkdata(i+ioff,j-1+joff);   u22 = blkdata(i+ioff,j+joff);   u23 = blkdata(i+ioff,j+1+joff);
                u31 = blkdata(i+1+ioff,j-1+joff); u32 = blkdata(i+1+ioff,j+joff); u33 = blkdata(i+1+ioff,j+1+joff);

                blkdata2(2*i-1,2*j-1) = (11.0*u11 + 88.0*u12 - 11.0*u13 - ...
                                            4.0*u21 - 32.0*u22 + 4.0*u23 + u31 + ...
                                            8.0*u32 - 1.0*u33) * 0.015625;
                blkdata2(2*i-1,2*j) = (11.0*u11 - 88.0*u12 - 11.0*u13 - ...
                                            4.0*u21 + 32.0*u22 + 4.0*u23 + u31 - ...
                                            8.0*u32 - 1.0*u33) * (-0.015625);
                blkdata2(2*i,2*j-1) = (5.0*u11 + 40.0*u12 - 5.0*u13 + ...
                                            4.0*u21 + 32.0*u22 - 4.0*u23 - 1.0*u31 - ...
                                            8.0*u32 + u33) * 0.015625;
                blkdata2(2*i,2*j) =  (5.0*u11 - 40.0*u12 - 5.0*u13 + ...
                                            4.0*u21 - 32.0*u22 - 4.0*u23 - 1.0*u31 + ...
                                            8.0*u32 + u33) * (-0.015625);

            end

            % upper left
            i = 1;
            j = nyblk;
            ioff = 1;
            joff = -1;

            u11 = blkdata(i-1+ioff,j-1+joff); u12 = blkdata(i-1+ioff,j+joff); u13 = blkdata(i-1+ioff,j+1+joff);
            u21 = blkdata(i+ioff,j-1+joff);   u22 = blkdata(i+ioff,j+joff);   u23 = blkdata(i+ioff,j+1+joff);
            u31 = blkdata(i+1+ioff,j-1+joff); u32 = blkdata(i+1+ioff,j+joff); u33 = blkdata(i+1+ioff,j+1+joff);

            blkdata2(2*i,2*j) = 4*((5*u11)/256 - (5*u12)/64 + (55*u13)/256 + u21/64 - u22/16 + ...
                            (11*u23)/64 - u31/256 + u32/64 - (11*u33)/256);
            blkdata2(2*i-1,2*j) = 4*((11*u11)/256 - (11*u12)/64 + (121*u13)/256 - u21/64 + ...
                            u22/16 - (11*u23)/64 + u31/256 - u32/64 + (11*u33)/256);
            blkdata2(2*i-1,2*j-1) = 4*((11*u12)/64 - (11*u11)/256 + (55*u13)/256 + u21/64 - u22/16 - ...
                            (5*u23)/64 - u31/256 + u32/64 + (5*u33)/256);
            blkdata2(2*i,2*j-1) = 4*((5*u12)/64 - (5*u11)/256 + (25*u13)/256 - u21/64 + u22/16 + ...
                            (5*u23)/64 + u31/256 - u32/64 - (5*u33)/256);

            % upper right
            i = nxblk;
            j = nyblk;
            ioff = -1;
            joff = -1;

            u11 = blkdata(i-1+ioff,j-1+joff); u12 = blkdata(i-1+ioff,j+joff); u13 = blkdata(i-1+ioff,j+1+joff);
            u21 = blkdata(i+ioff,j-1+joff);   u22 = blkdata(i+ioff,j+joff);   u23 = blkdata(i+ioff,j+1+joff);
            u31 = blkdata(i+1+ioff,j-1+joff); u32 = blkdata(i+1+ioff,j+joff); u33 = blkdata(i+1+ioff,j+1+joff);

            blkdata2(2*i,2*j) = (u11 - 4.0*u12 + 11.0*u13 - 4.0*u21 + ...
                                        16.0*u22 - 44.0*u23 + 11.0*u31 - 44.0*u32 ...
                                        + 121.0*u33) * 0.015625;
            blkdata2(2*i-1,2*j) = (u11 - 4.0*u12 + 11.0*u13 - 4.0*u21 + ...
                                        16.0*u22 - 44.0*u23 - 5.0*u31 + 20.0*u32 - ...
                                        55.0*u33) * (-0.015625);
            blkdata2(2*i-1,2*j-1) = (u11 - 4.0*u12 - 5.0*u13 - 4.0*u21 + ...
                                        16.0*u22 + 20.0*u23 - 5.0*u31 + 20.0*u32 + ...
                                        25.0*u33) * 0.015625;
            blkdata2(2*i,2*j-1) = (u11 - 4.0*u12 - 5.0*u13 - 4.0*u21 + ...
                                        16.0*u22 + 20.0*u23 + 11.0*u31 - 44.0*u32 ...
                                        - 55.0*u33) * (-0.015625);

            % lower right
            i = nxblk;
            j = 1;
            ioff = -1;
            joff = 1;

            u11 = blkdata(i-1+ioff,j-1+joff); u12 = blkdata(i-1+ioff,j+joff); u13 = blkdata(i-1+ioff,j+1+joff);
            u21 = blkdata(i+ioff,j-1+joff);   u22 = blkdata(i+ioff,j+joff);   u23 = blkdata(i+ioff,j+1+joff);
            u31 = blkdata(i+1+ioff,j-1+joff); u32 = blkdata(i+1+ioff,j+joff); u33 = blkdata(i+1+ioff,j+1+joff);

            blkdata2(2*i-1,2*j-1) = (11.0*u11 - 4.0*u12 + u13 - 44.0*u21 + ...
                            16.0*u22 - 4.0*u23 - 55.0*u31 + 20.0*u32 - ...
                            5.0*u33) * (-0.015625);
            blkdata2(2*i-1,2*j) = (5.0*u11 + 4.0*u12 - 1.0*u13 - ...
                            20.0*u21 - 16.0*u22 + 4.0*u23 - 25.0*u31 - ...
                            20.0*u32 + 5.0*u33) * (-0.015625);
            blkdata2(2*i,2*j-1) = (11.0*u11 - 4.0*u12 + u13 - 44.0*u21 + ...
                            16.0*u22 - 4.0*u23 + 121.0*u31 - 44.0*u32 ...
                            + 11.0*u33) * (0.015625);
            blkdata2(2*i,2*j) = (5.0*u11 + 4.0*u12 - 1.0*u13 - ...
                            20.0*u21 - 16.0*u22 + 4.0*u23 + 55.0*u31 + ...
                            44.0*u32 - 11.0*u33) * (0.015625);

            % lower left
            i = 1;
            j = 1;
            ioff = 1;
            joff = 1;

            u11 = blkdata(i-1+ioff,j-1+joff); u12 = blkdata(i-1+ioff,j+joff); u13 = blkdata(i-1+ioff,j+1+joff);
            u21 = blkdata(i+ioff,j-1+joff);   u22 = blkdata(i+ioff,j+joff);   u23 = blkdata(i+ioff,j+1+joff);
            u31 = blkdata(i+1+ioff,j-1+joff); u32 = blkdata(i+1+ioff,j+joff); u33 = blkdata(i+1+ioff,j+1+joff);

            blkdata2(2*i-1,2*j-1) = (11.0*u13 - ...
                            44.0*u21 + 16.0*u22 - 4.0*u23 + 11.0*u31 - ...
                            4.0*u32 + u33 - 44.0*u12 + 121.0*u11) * 0.015625;
            blkdata2(2*i-1,2*j) = (55.0*u11 + 44.0*u12 - 11.0*u13 - ...
                            20.0*u21 - 16.0*u22 + 4.0*u23 + 5.0*u31 + ...
                            4.0*u32 - 1.0*u33) * 0.015625;
            blkdata2(2*i,2*j-1) = (55.0*u11 - 20.0*u12 + 5.0*u13 + ...
                            44.0*u21 - 16.0*u22 + 4.0*u23 - 11.0*u31 + ...
                            4.0*u32 - 1.0*u33) * 0.015625;
            blkdata2(2*i,2*j) = (- 5.0*u13 + ...
                            20.0*u21 + 16.0*u22 - 4.0*u23 - 5.0*u31 - ...
                            4.0*u32 + u33 + 25.0*u11 + 20.0*u12) * 0.015625;

        end

        % copy variable
        blkdata = blkdata2;

      end

      % get size of this newly created block
      [nxblk nyblk] = size(blkdata);

      % determine the starting (i,j) indices in global array
      iglb = int32((xlo - xmin) / dxmin + 1);
      jglb = int32((ylo - ymin) / dymin + 1);

      % copy data
      for i = 0:nxblk-1
        for j = 0:nyblk-1
          vdata(jglb+j,iglb+i) = blkdata(i+1,j+1);
        end
      end

    end

  end

end
