%% matrix.m --- local nonlinear transient analysis
% |2011-12-07, wyyjtu@gmail.com|
% |2012-12-25, add nonlinear iteration.|
%
% MATLAB R2012a
% ANSYS 12.1
%
% _code:_

%% matrix
% The process for a typical transient analysis.
%
% 1. Building the Model.  2. Applying Loads and Obtaining the
% Solution.  3. Reviewing the Results.
function matrix(varargin)
    nvarargs = length(varargin);
    if nvarargs == 3 && strcmp(varargin{1}, 'post')
        node = varargin{2};
        dof = varargin{3};
        load('result.mat')
        if ~ismember(node, map_node) || ~ismember(dof, map_dof)
            fprintf('NODE or DOF error\n')
            return
        end
        resplot(node, dof, tr, map_node, map_dof, rst, time_s)
        return
    elseif nvarargs <= 2
        opt = {'linear', 'nonlinear', 'irs', 'lr'};
        for i = 1:nvarargs
            if ~ismember(varargin{i}, opt)
                type('README.md')
                return
            end
        end
        s = ismember(opt, varargin);
        opt_linear = s(1);
        opt_irs = s(3);
        if s(2) == 0
            opt_linear = 1;
        end
        if s(4) == 0
            opt_irs = 1;
        end
    else
        type('README.md')
        return
    end

    % ANSYS
    global ansys
    ansys = '"C:\Program Files\ANSYS Inc\v121\ansys\bin\WINX64\ansys121" -b';

    % Model
    system([ansys ' -i model.mac -o vm1.out']);
    fprintf('READ %s\n', 'K.TXT', 'M.TXT', 'K.mapping');
    [k, rhs] = readhb('K.TXT');
    [m, r0] = readhb('M.TXT');
    [map_node, map_dof] = readmap('K.mapping', size(k,1));
    % Model reduction
    fprintf('\nNatural frequencies (Hz)\n');
    freq(k, m, 'EXACT');
    if opt_irs == 1
        [mr, kr, tr] = irs(k, m, map_node);
    else
        [mr, kr, tr] = lrbm(k, m, map_node, map_dof);
    end

    tic
    % Apply loads.                      file1 file2
    % Specify load options [alpha, beta, dof1, dof2, scale]
    fprintf('\nApply loads and solve\n');
    [cr, fs, time_s, step] = loads(mr, kr, tr, map_dof, 'I-ELC180.AT2', ...
                                  'I-ELC270.AT2', [0.6, 0.006, 3, 1, 1]);
    % Solve
    if opt_linear == 1
        rst = hpdl(mr, cr, kr, tr, fs, step, map_node, map_dof, 'linear');
    else
        rst = hpdl(mr, cr, kr, tr, fs, step, map_node, map_dof, 'nonlinear');
    end
    toc
    rst = rst(1:size(kr,1), :);
    save('result.mat', 'tr', 'map_node', 'map_dof', 'rst', 'time_s')
end

%% readhb
% Read matrix in Harwell-Boeing format.
function [matrix, rhs] = readhb(file_name)
    fid = fopen(file_name);
    % HB file header
    tline = fgetl(fid);
    tline = fgetl(fid);
    b = sscanf(tline, '%d', 5);
    tline = fgetl(fid);
    c = sscanf(tline, '%*s%d%d%d%d');
    tline = fgetl(fid);
    if b(5) ~= 0
        tline = fgetl(fid);
    end
    % matrix data
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
% Read the mapping file.
function [map_node, map_dof] = readmap(file_name, n)
    map_node = zeros(n,1);
    map_dof = zeros(n,1);
    fid = fopen(file_name);
    tline = fgetl(fid);
    for i = 1:n
        map_node(i,1) = fscanf(fid, '%*d%d');
        map_dof(i,1) = u2dof(fscanf(fid, '%s', 1));
    end
    fclose(fid);
end

%% u2dof
% U to DOF.
function dof = u2dof(u)
    switch u
      case 'UX'
        dof = 1;
      case 'UY'
        dof = 2;
      case 'UZ'
        dof = 3;
      case 'ROTX'
        dof = 4;
      case 'ROTY'
        dof = 5;
      case 'ROTZ'
        dof = 6;
      otherwise
        dof = 0;
        fprintf('\n*** u2dof ERROR ***\n\n');
    end
end

%% freq
% Natural frequencies.
function frequencies = freq(k, m, flags)
    n = 15;
    d = eigs(k, m, n, 'sm');
    d = sort(d, 1);
    % units: Hz
    frequencies = sqrt(d) / (2 * pi);
    fprintf([sprintf('%7s:', flags), sprintf('%11.5f', frequencies), '\n']);
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
function [mr, kr, t] = irs(k, m, map_node)
    [meqid, seqid, nmdof] = mseqid(map_node);
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
        % frequencies, error
        freqs = freq(kr, mr, sprintf('IRS-%d',i));
        if i > 0
            error = max(abs((freqs - freqs_old) ./ freqs_old));
            if error < 0.001
                break
            end
        end
        freqs_old = freqs;
        % the transformation matrix for IRS
        i = i + 1;
        temp = msm + mss * t(seqid, :) * (mr \ kr);
        t(seqid, :) = t_static + (kss \ temp);
    end
end

%% mseqid
% Get meqid, seqid, nmdof.
function [meqid, seqid, nmdof] = mseqid(map_node)
    nd_list = load('MNODES.TXT');
    meqid = ismember(map_node, nd_list);
    seqid = not(meqid);
    nmdof = size(find(meqid), 1);
end

%% lrbm
% Local rigid body mode.
%
% The transformation matrix is
%
% $$\mathbf{T}=\mathbf{K}^{-1}(\mathbf{R}^T)^+\mathbf{R}^T\mathbf{KR}.$$
function [mr, kr, t] = lrbm(k, m, map_node, map_dof)
    r = getr(map_node, map_dof, size(k,1));
    % t
    prkr = pinv(r') * (r' * k * r);
    t = k \ prkr;
    % reduced K, M
    mr = t' * m * t;
    kr = t' * k * t;
    % frequencies
    freq(kr, mr, 'LR');
end

%% getr
% Computer r.
function r = getr(map_node, map_dof, n)
    nreg = load('RIGID.TXT');
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
            % fill r
            indexr = map_node == nd_list(j);
            indext = map_dof(indexr, 1);
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
function rst = hpdl(m, c, k, trans, f, step, map_node, map_dof, flags)
    n = size(k, 1);
    % H
    im = inv(m);
    h = zeros(2*n, 2*n);
    h(1:n, n+1:end) = eye(n);
    h(n+1:end, 1:n) = - im * k;
    h(n+1:end, n+1:end) = - im * c;
    ih = inv(h);
    % T
    dt = step / 2^20;
    ta = h*dt + (h*dt)^2/2 + (h*dt)^3/6 + (h*dt)^4/24;
    for i = 1:20
        ta = 2*ta + ta*ta;
    end
    t = eye(size(h)) + ta;
    % f
    f = im * f;
    % hpdl iteration
    if strcmp(flags, 'linear')
        rst = hpdl_iteration(ih, t, f, step);
    elseif strcmp(flags, 'nonlinear')
        rst = hpdl_non_iteration(ih, t, f, step, im, trans, map_node);
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
function rst = hpdl_non_iteration(ih, t, f, step, im, trans, map_node)
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
    non_nodes = load('NONL.TXT');
    non_d = zeros(size(non_nodes,1)*3, 1);
    non_iter = 0;
    ff = zeros(size(trans, 1), 1);
    non_print = '[%6d/%6d] %2d, Cumulative Iteration Number: %7d\n';
    non_p_len = length(sprintf(non_print, 0, 0, 0, 0));
    non_pbs = char(reshape([92;98]*ones(1,non_p_len), 1, []));
    % initial nonlinear analysis
    fprintf(non_print, 0, 0, 0, 0);
    system([ansys ' -i nonlinear.mac -o vm1.out']);
    for i = 1:nstep
        [status, result] = system('move file.r002 file.r001');
        for j = 0:99
            non_f_trans = trans' * ff;
            r(n+1:end, 1) = f(:, i);
            r1(n+1:end, 1) = (f(:, i+1) - im * non_f_trans - f(:, i)) / step;
            v_temp = t*(v + ih*(r + ih*r1)) - ih*(r + ih*r1 + r1*step);
            % error
            if j > 0
                error = max(abs((v_temp - v_old) ./ v_old));
                if j == 20
                    str = [non_pbs, sprintf('Nonconvergence: step %d\n', i)];
                    break
                elseif error < 0.001
                    str = non_pbs;
                    break
                end
            end
            v_old = v_temp;
            % output degree-of-freedom constraints at nodes
            v_trans = trans * v_temp(1:n, 1);
            for ii = 1:size(non_nodes,1)
                index = find(map_node == non_nodes(ii));
                non_d(ii*3-2:ii*3, 1) = v_trans(index(1:3), 1);
            end
            fid = fopen('D.TXT', 'w');
            for ii = 1:size(non_d,1)
                fprintf(fid, '%25.15E\n', non_d(ii,1));
            end
            fclose(fid);
            % restart a previous analysis
            fid = fopen('RESTART.TXT', 'w');
            fprintf(fid, ['FINISH$/CLEAR$/FILNAME,file\n' ...
                          '/SOLU$ANTYPE,,REST,%d\n' ...
                          'TIME,%g$nonrest\n'], i, i*step);
            fclose(fid);
            system([ansys ' -i RESTART.TXT -o vm1.out']);
            % load nodal reaction forces
            non_f = load('F.TXT');
            ff = zeros(size(trans, 1), 1);
            for ii = 1:size(non_nodes,1)
                index = find(map_node == non_nodes(ii));
                ff(index(1:3), 1) = ff(index(1:3), 1) + non_f(ii*3-2:ii*3, 1);
            end
        end
        non_iter = non_iter + j;
        fprintf([str, non_print], i, nstep, j, non_iter);
        % result
        v = v_temp;
        f(:, i+1) = f(:, i+1) - im * non_f_trans;
        rst(:, i) = v;
    end
end

%% loads
% Apply acceleration loads in a transient analysis.
function [c, fs, time_s, step] = loads(m, k, t, map_dof, file1, file2, opt)
    alpha = opt(1);
    beta = opt(2);
    dof1 = opt(3);
    dof2 = opt(4);
    scale = opt(5);
    % damping
    c = alpha * m + beta * k;
    % load earthquake
    [step, acce_1] = load_earthquake(file1);
    [step, acce_2] = load_earthquake(file2);
    % direction
    index_1 = map_dof == dof1;
    index_2 = map_dof == dof2;
    % transformed force sequence
    fs_1 = - m * pinv(t) * index_1 * acce_1' * 9.81;
    fs_2 = - m * pinv(t) * index_2 * acce_2' * 9.81;
    fs = fs_1 + fs_2;
    % scaling
    step = step * scale;
    fs = fs(:, scale:scale:end);
    % time_s
    time = step * size(fs, 2);
    time_s = step:step:time;
end

%% load_earthquake
% Load ground motion data.
%
% The Pacific Earthquake Engineering Research Center (PEER)
% ground motion database: http://peer.berkeley.edu/smcat/
function [step, acce] = load_earthquake(file_name)
    fid = fopen(file_name);
    for i = 1:4
        tline = fgetl(fid);
    end
    step = sscanf(tline, 'NPTS= %*d, DT= %f SEC');
    acce = fscanf(fid, '%f');
    fclose(fid);
end

%% resplot
% Plot and save nodal DOF result.
function resplot(node, dof, t, map_node, map_dof, rst, time_s)
    index = map_node == node & map_dof == dof;
    t_nd = t(index, :);
    y = t_nd * rst;
    plot(time_s, y)
    % save
    outmat = [time_s; y]';
    save(sprintf('rst-%d-%d.txt', node, dof), 'outmat', '-ASCII')
end

% matrix.m ends here
