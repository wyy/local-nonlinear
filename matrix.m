%% matrix.m --- local nonlinear transient analysis
% |2011-12-07, Wang Yingying|
%
% |modified: add latex (2012-06-28, wyy).|
% |modified: earthquake analysis (2012-08-05, wyy).|
%
% *Matlab version (R2011, R2012)*.
%
% _code:_

%% matrix
% Verify the *path to ANSYS*.
%
% Customize *target error* of IRS method.
function matrix()
    clc
    % settings
    global ansys
    ansys = '"C:\Program Files\ANSYS Inc\v121\ansys\bin\WINX64\ansys121" -b';
    % call ANSYS
    fprintf('call ANSYS\n\n');
    [status, result] = system('clean');
    system([ansys ' -i model.mac -o vm1.out']);

    % read [k] [m] and [mapping]
    [k, rhs] = readhb('K_RHS.TXT');
    [m, r0] = readhb('M.TXT');
    [dofmp, dofmpxyn] = readmap('K_RHS.mapping', size(k,1));
    % plots the sparsity pattern of k
    % spy(k, 5)
    % exact solution
    fprintf('\nEXACT SOLUTION\n');
    eigs_vector = get_eigs(k, m, 'EXACT');
    % model reduction
    [ml, kl, tl] = lrbm(k, m, dofmp, dofmpxyn);
    [mi, ki, ti] = irs(k, m, dofmp);
    % transient analysis
    tic
    fprintf('\nTRANSIENT ANALYSIS WITH HPD-L\n');
    trans_hpd(mi, ki, ti, rhs, dofmp, dofmpxyn, 'irs')
    % trans_hpd(ml, kl, tl, rhs, dofmp, dofmpxyn, 'lr')
    % trans_hpd(mi, ki, ti, rhs, dofmp, dofmpxyn, 'resume')
    toc
end

%% get_eigs
% Get the n smallest magnitude eigenvalues, which is displayed on the
% screen and saved to `log-eigs.txt' at the same time.
%
% Customize *n*.
function eigs_vector = get_eigs(k, m, flags)
    n = 12;
    d = eigs(k, m, n, 'sm');
    d = sort(d, 1);
    eigs_vector = sqrt(d) / (2*pi);
    % write to log file
    fid = fopen('log-eigs.txt', 'a');
    fprintf(fid, format_out_str(flags, eigs_vector));
    fclose(fid);
end

%% format_out_str
% Get formated out_string, and print it on the screen.
%
% Customize *display format*.
function out_string = format_out_str(flags, eigs_v)
    out_string = [sprintf('%6s:', flags), sprintf('%11.5f', eigs_v), '\n'];
    fprintf(out_string);
end

%% readhb
% Read matrix from file with Harwell-Boeing format.
function [matrix, rhs] = readhb(file_name)
    fprintf('READ %s\n', file_name);
    fid = fopen(file_name);
    % read header from HB file
    tline = fgetl(fid);
    tline = fgetl(fid);
    b = sscanf(tline, '%d', 5);
    tline = fgetl(fid);
    c = sscanf(tline, '%*s%d%d%d%d');
    tline = fgetl(fid);
    if b(5) ~= 0
        tline = fgetl(fid);
        % d = sscanf(tline, '%*s%d%d');
    end
    % read matrix data
    ncolu = b(2); % number of line(column) info
    nrow = b(3); % number of row info
    nrhs = b(5); % number of RHS
    nr = c(1); % number of rows
    nc = c(2); % number of columns
    columninfo = fscanf(fid, '%d', ncolu); % column info
    rowindex = fscanf(fid, '%d', nrow); % row index
    temp = fscanf(fid, '%gD%d', [2, nrow])'; % matrix elements
    matrixele = temp(:,1) .* 10.^temp(:,2);
    temp = fscanf(fid, '%gD%d', [2, nrhs])'; % RHS
    rhs = temp(:,1) .* 10.^temp(:,2);
    fclose(fid);
    % assemble matrix
    colindex = zeros(nrow, 1);
    for i = 1:(ncolu-1)
        colindex(columninfo(i) : (columninfo(i+1)-1)) = i;
    end
    colidx = [colindex; rowindex];
    rowidx = [rowindex; colindex];
    matele = [matrixele; matrixele];
    matrix = sparse(rowidx, colidx, matele, nr, nc);
    for i = 1:nr
        matrix(i,i) = matrix(i,i) / 2;
    end
end

%% readmap
% Read filename.mapping file.
function [dofmp, dofmpxyn] = readmap(file_name, n)
    fprintf('READ %s\n', file_name);
    dofmp = zeros(n,1);
    dofmpxyn = zeros(n,1);
    fid = fopen(file_name);
    tline = fgetl(fid);
    for i = 1:n
        dofmp(i,1) = fscanf(fid, '%*d%d');
        dofmpxyn(i,1) = uxy2num(fscanf(fid, '%s', 1));
    end
    fclose(fid);
end

%% uxy2num
% uxy to number.
function num = uxy2num(uxy)
    switch uxy
      case 'UX'
        num = 1;
      case 'UY'
        num = 2;
      case 'UZ'
        num = 3;
      case 'ROTX'
        num = 4;
      case 'ROTY'
        num = 5;
      case 'ROTZ'
        num = 6;
      otherwise
        num = 0;
        fprintf('\n*** UXY2NUM error ***\n\n');
    end
end

%% irs
% IRS method.
%
% The equation of motion of the structure
%
% $$\left[\!\! \begin{array}{cc}
%     \mathbf{M}_{mm} & \mathbf{M}_{ms} \\
%     \mathbf{M}_{sm} & \mathbf{M}_{ss}
%   \end{array} \!\!\right] \left\{\!\!\!\! \begin{array}{c}
%   \ddot{\mathbf{x}}_m \\\ddot{\mathbf{x}}_s
% \end{array} \!\!\!\! \right \} + \left[ \!\!\begin{array}{cc}
%   \mathbf{K}_{mm} & \mathbf{K}_{ms} \\
%   \mathbf{K}_{sm} & \mathbf{K}_{ss}
%   \end{array} \!\!\right] \left\{\!\!\!\! \begin{array}{c}
%   \mathbf{x}_m \\\mathbf{x}_s
% \end{array} \!\!\!\! \right \} = \left \{\!\!\!\! \begin{array}{c}
%   \mathbf{f}_{m} \\ \mathbf{0}
% \end{array} \!\!\!\! \right \}.$$
%
% _Static reduction_ :
% Eliminate the slave degrees of freedom
%
% $$\left\{\!\!\!\! \begin{array}{c}
%    \mathbf{x}_m \\\mathbf{x}_s
%  \end{array} \!\!\!\! \right \} = \left[\!\! \begin{array}{c}
%    \mathbf{I} \\
%    -\mathbf{K}_{ss}^{-1} \mathbf{K}_{sm}
%  \end{array} \!\!\right]\{\mathbf{x}_m\} = \mathbf{T}_s\mathbf{x}_m.$$
%
% _Iterated IRS techniques_ :
% The calculation involves only the lower part of the transformation
%
% $$\mathbf{T}_{IRS,i} = \left[\begin{array}{c}
%    \mathbf{I} \\ \mathbf{t}_{IRS,i} \end{array}\right],$$
%
% where
%
% $$\mathbf{t}_{IRS,i+1} = -\mathbf{K}_{ss}^{-1}\mathbf{K}_{sm} +
% \mathbf{K}_{ss}^{-1}(\mathbf{M}_{sm}+\mathbf{M}_{ss}\mathbf{t}_{IRS,i})
% \mathbf{M}_{IRS,i}^{-1}\mathbf{K}_{IRS,i}.$$
function [mr, kr, t] = irs(k, m, dofmp)
    fprintf('\nBEGIN IRS METHOD\n');
    % split K, M into sub matrices relating to the master degrees of
    % freedom
    [meqid, seqid, nmdof] = mseqid(dofmp);
    kss = k(seqid, seqid);
    ksm = k(seqid, meqid);
    mss = m(seqid, seqid);
    msm = m(seqid, meqid);
    % the transformation matrix
    t = zeros(size(k, 1), nmdof);
    t(meqid, :) = eye(nmdof);
    % static reduction (Guyan)
    t_static = - (kss \ ksm);
    t(seqid, :) = t_static;
    % iterated IRS techniques
    i = 0;
    while 1
        mr = t' * m * t;
        kr = t' * k * t;
        % eigs, error
        eigs_vector = get_eigs(kr, mr, sprintf('IRS%3d',i));
        if i > 0
            error = max(abs((eigs_vector - eigs_old) ./ eigs_old));
            if error < 0.001
                break
            end
        end
        eigs_old = eigs_vector;
        % the transformation matrix for IRS
        i = i + 1;
        temp = msm + mss * t(seqid, :) * (mr \ kr);
        t(seqid, :) = t_static + (kss \ temp);
    end
end

%% mseqid
% Get meqid, seqid, nmdof, restor.
function [meqid, seqid, nmdof] = mseqid(dofmp)
    mnd_list = load('M_NODES.TXT');
    meqid = ismember(dofmp, mnd_list);
    seqid = not(meqid);
    nmdof = size(find(meqid), 1);
end

%% lrbm
% Local rigid body mode.
%
% The transformation matrix is
%
% $$\mathbf{T}=\mathbf{K}^{-1}(\mathbf{R}^T)^+\mathbf{R}^T\mathbf{KR}.$$
function [mr, kr, t] = lrbm(k, m, dofmp, dofmpxyn)
    fprintf('\nBEGIN LOCAL RIGID BODY MODE\n');
    r = getr(dofmp, dofmpxyn, size(k,1));
    % t
    prkr = pinv(r') * (r' * k * r);
    t = k \ prkr;
    % reduced K, M
    mr = t' * m * t;
    kr = t' * k * t;
    % eigs
    eigs_vector = get_eigs(kr, mr, 'LR');
end

%% getr
% Computer r.
function r = getr(dofmp, dofmpxyn, n)
    nreg = load('AMOUNT.TXT');
    r = zeros(n, 6*nreg);
    temp = eye(6);
    for i = 1:nreg
        nd_list = load(['NODES_' num2str(i) '.TXT']);
        ndxy = load(['NODESXY_' num2str(i) '.TXT']);
        ndxy = reshape(ndxy, size(nd_list,1), 3);
        xym = sum(ndxy) / size(ndxy,1);
        xm = xym(1); ym = xym(2); zm = xym(3);
        for j = 1:size(nd_list,1)
            x = ndxy(j,1); y = ndxy(j,2); z = ndxy(j,3);
            xd = x - xm; yd = y - ym;zd = z - zm;
            temp(1,5) = zd;
            temp(1,6) = - yd;
            temp(2,4) = - zd;
            temp(2,6) = xd;
            temp(3,4) = yd;
            temp(3,5) = - xd;
            % fill to r
            indexr = dofmp == nd_list(j);
            indext = dofmpxyn(indexr, 1);
            r(indexr, (6*(i-1)+1):(6*i)) = temp(indext, :);
        end
    end
end

%% hpdl
% High Precision Direct Integration - Linear.
%
% Convert
%
% $$[\mathbf{M}]\{\ddot{\mu}\} + [\mathbf{C}]\{\dot{\mu}\} +
% [\mathbf{K}]\{\mu\} = \{f(t)\}$$
%
% to
%
% $$\{\dot{\nu}\}=[\mathbf{H}]\{\nu\}+\{r(t)\},$$
%
% where
%
% $$\{\nu\}^T=\{\mu^T,\dot{\mu}^T\}.$$
%
% Dimensions of matrices
%
% $$m,k,c \in \mathbf{R}^{n \times n}, \quad f \in \mathbf{R}^{n \times
%   (nstep+1)}, \quad v \in \mathbf{R}^{2n \times 1}, \quad rst \in
% \mathbf{R}^{2n \times nstep}.$$
function rst = hpdl(m, c, k, f, step, trans, dofmp, dofmpxyn, varargin)
    n = size(k, 1);
    % H
    im = inv(m);
    h = zeros(2*n, 2*n);
    h(1:n, n+1:end) = eye(n);
    h(n+1:end, 1:n) = - im * k;
    h(n+1:end, n+1:end) = - im * c;
    ih = inv(h);
    % computer T
    dt = step / 2^20;
    ta = h*dt + (h*dt)^2/2 + (h*dt)^3/6 + (h*dt)^4/24;
    for i = 1:20
        ta = 2*ta + ta*ta;
    end
    t = eye(size(h)) + ta;
    % f
    f = im * f;
    % hpdl iteration
    if length(varargin) == 0
        rst = hpdl_iteration(ih, t, f, step);
    else
        % iteration with nonlinear force
        rst = hpdl_non_iteration(im, ih, t, f, step, trans, dofmp, dofmpxyn);
    end
end

%% hpdl_iteration
% hpdl iteration
function rst = hpdl_iteration(ih, t, f, step)
    n = size(f, 1);
    nstep = size(f, 2);
    r = zeros(2*n, 1);
    r1 = zeros(2*n, 1);
    rst = zeros(2*n, nstep);
    % initilize
    f = [zeros(n,1), f];
    v = zeros(2*n, 1);
    for i = 1:nstep
        r(n+1:end, 1) = f(:, i);
        r1(n+1:end, 1) = (f(:, i+1) - f(:, i)) / step;
        v = t*(v + ih*(r + ih*r1)) - ih*(r + ih*r1 + r1*step);
        rst(:, i) = v;
    end
end

%% hpdl_non_iteration
% hpdl iteration with nonlinear force
function rst = hpdl_non_iteration(im, ih, t, f, step, trans, dofmp, dofmpxyn)
    n = size(f, 1);
    nstep = size(f, 2);
    r = zeros(2*n, 1);
    r1 = zeros(2*n, 1);
    rst = zeros(2*n, nstep);
    % initilize
    f = [zeros(n,1), f];
    v = zeros(2*n, 1);
    % nonlinear related variables
    global ansys
    non_f_full = zeros(size(trans, 1), 1);
    non_nodes = load('NONL.DAT');
    non_d = zeros(size(non_nodes,1)*3, 1);
    non_f = zeros(size(non_d));
    fid = fopen('NON.TXT', 'w');
    fprintf(fid, '%10d\n', size(non_nodes,1));
    fclose(fid);
    system([ansys ' -i nonlinear.mac -o vm1.out']);
    for i = 1:nstep
        fprintf('step %6d of %d, ', i, nstep);
        [status, result] = system('move FILE.r002 FILE.r001');
        for j = 0:99
            non_f_trans = trans' * non_f_full;
            r(n+1:end, 1) = f(:, i);
            r1(n+1:end, 1) = (f(:, i+1) - im * non_f_trans - f(:, i)) / step;
            v_temp = t*(v + ih*(r + ih*r1)) - ih*(r + ih*r1 + r1*step);
            % error
            if j > 0
                error = max(abs((v_temp - v_old) ./ v_old));
                if j == 20
                    fprintf('nonconvergence, ');
                    break
                elseif error < 0.001
                    break
                end
            end
            v_old = v_temp;
            % nonlinear force
            % output degree-of-freedom constraints at nodes
            v_trans = trans * v_temp(1:n, 1);
            for ii = 1:size(non_nodes,1)
                index = find(dofmp == non_nodes(ii));
                non_d(ii*3-2:ii*3, 1) = v_trans(index(1:3), 1);
            end
            fid = fopen('D.TXT', 'w');
            for ii = 1:size(non_d,1)
                fprintf(fid, '%25.15E\n', non_d(ii,1));
            end
            fclose(fid);
            % restart a previous analysis
            fid = fopen('RESTART.TXT', 'w');
            fprintf(fid, ['FINISH$/CLEAR$/FILNAME,FILE\n' ...
                          '/SOLU$ANTYPE,,REST,%d\n' ...
                          'TIME,%g$nonrest\n'], i, i*step);
            fclose(fid);
            system([ansys ' -i RESTART.TXT -o vm1.out']);
            % load nodal reaction forces
            non_f = load('F.TXT');
            for ii = 1:size(non_nodes,1)
                index = find(dofmp == non_nodes(ii));
                non_f_full(index(1:3), 1) = non_f(ii*3-2:ii*3, 1);
            end
        end
        v = v_temp;
        f(:, i+1) = f(:, i+1) - im * non_f_trans;
        rst(:, i) = v;
        fprintf('nonlinear iteration: %d\n', j);
        fid = fopen('log-iteration.txt', 'a');
        fprintf(fid, '%d\n', j);
        fclose(fid);
    end
end

%% trans_hpd
% Transient analysis with HPD-L.
function trans_hpd(m, k, t, rhs, dofmp, dofmpxyn, method)
    n = size(k, 1);
    % damping
    alpha = 5;
    beta = 0;
    c = alpha * m + beta * k;
    % load earthquake
    [step, asequence_x] = rd_earthquake('I-ELC180.AT2');
    [step, asequence_y] = rd_earthquake('I-ELC270.AT2');
    % direction x y
    index_x = dofmpxyn == 1;
    index_y = dofmpxyn == 2;
    % transformed force sequence
    fs_x = - m * pinv(t) * index_x * asequence_x' * 9.81;
    fs_y = - m * pinv(t) * index_y * asequence_y' * 9.81;
    fs = fs_x + fs_y;
    % modify step
    scale = 2;
    step = step * scale;
    fs = fs(:, scale:scale:end);

    if strcmp(method, 'resume')
        % resume
        rst = load(['log-rst.txt']);
    else
        % hpdl
        rst = hpdl(m, c, k, fs, step, t, dofmp, dofmpxyn, 'non');
        save(['log-rst.txt'], 'rst', '-ASCII')
    end
    % plot displacement
    time = step * size(fs, 2);
    times = step:step:time;
    resplot(1546, 1, full(t), dofmp, dofmpxyn, times, rst(1:n,:), method)
end

%% rd_earthquake
% Read acceleration time history and time step.
%
% The Pacific Earthquake Engineering Research Center (PEER)
% ground motion database: http://peer.berkeley.edu/smcat/
function [step, asequence] = rd_earthquake(file_name)
    fid = fopen(file_name);
    for i = 1:4
        tline = fgetl(fid);
    end
    step = sscanf(tline, 'NPTS= %*d, DT= %f SEC');
    asequence = fscanf(fid, '%f');
    fclose(fid);
    save(['log-' file_name '.txt'], 'asequence', '-ASCII')
end

%% resplot
% Plot and save the result of a node (node_num: ndof).
function resplot(node_num, ndof, t, dofmp, dofmpxyn, times, rst, method)
    index = dofmp == node_num & dofmpxyn == ndof;
    t_res = t(index, :);
    y = t_res * rst;
    plot(times, y)
    hold on
    % save in files
    times = times';
    y = y';
    save('log-times.txt', 'times', '-ASCII')
    save(['log-' method '-nsolu.txt'], 'y', '-ASCII')
end

% matrix.m ends here
