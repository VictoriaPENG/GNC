function topo = generate_topology(cfg)
% ============================================================
% GENERATE_TOPOLOGY
% 多模式传感器网络拓扑生成函数
%
% 本函数用于生成森林防火无线传感器网络的节点空间拓扑，
% 支持多种典型部署方式，用于后续多跳路由与通信仿真。
%
% ---------------- 支持的拓扑模式 ----------------
% 1) uniform : 随机均匀布点
% 2) cluster : 分簇布点（多个局部高密度区域）
% 3) road    : 沿道路/走廊布点（带横向扩散）
%
% 同时支持"一个区域很多传感器"（热点区域叠加）
%
% ---------------- 输入参数 cfg ----------------
% 必选字段：
%   cfg.N          : 传感器节点数量（不含基站）
%   cfg.Lx, cfg.Ly: 监测区域尺寸
%   cfg.topo_mode : 'uniform' | 'cluster' | 'road'
%
% 常用可选字段：
%   cfg.seed       : 随机种子（保证可复现）
%   cfg.BS_pos     : 基站位置 [x y]
%
% 热点（高密度区域）参数：
%   cfg.hotspot_enable : 是否启用热点
%   cfg.hotspot_ratio  : 热点节点占比
%   cfg.hotspot_sigma  : 热点内扩散尺度
%   cfg.hotspot_center : 热点中心（默认火区）
%
% 分簇参数：
%   cfg.Kc                         : 簇数量
%   cfg.cluster_sigma              : 簇内节点扩散尺度
%   cfg.cluster_center_mode        : 簇中心生成方式
%   cfg.cluster_mix_uniform_ratio  : 混入均匀节点比例
%
% 道路参数：
%   cfg.road_count                 : 道路数量
%   cfg.road_width                 : 道路宽度（横向高斯扩散）
%   cfg.road_mix_uniform_ratio     : 混入均匀节点比例
%
% ---------------- 输出 topo ----------------
% topo.pos_sensors : N x 2 传感器节点坐标
% topo.pos_all     : (N+1) x 2 所有节点坐标（最后一行为 BS）
% topo.BS_pos      : 基站坐标
% topo.meta        : 拓扑元信息（用于可视化/调试）
%
% ============================================================

    % ----------- 参数完整性检查 -----------
    mustHave(cfg, {'N','Lx','Ly','topo_mode'});

    N  = cfg.N;
    Lx = cfg.Lx;
    Ly = cfg.Ly;

    % ----------- 固定随机种子，保证可复现 -----------
    if isfield(cfg,'seed') && ~isempty(cfg.seed)
        rng(cfg.seed,'twister');
    end

    % ----------- 基站位置 -----------
    if isfield(cfg,'BS_pos') && ~isempty(cfg.BS_pos)
        BS_pos = cfg.BS_pos(:)';
    else
        BS_pos = [Lx/2, Ly/2];   % 默认放在区域中心
    end

    topo_mode = lower(string(cfg.topo_mode));

    % 用于记录拓扑附加信息（画图/调试/论文说明）
    meta = struct();

    % 传感器节点位置（N x 2）
    pos = zeros(N,2);

    % =========================================================
    %               不同拓扑模式生成
    % =========================================================
    switch topo_mode

        % ---------------- 随机均匀布点 ----------------
        case "uniform"
            pos = uniform_points(N, Lx, Ly);
            meta.mode = 'uniform';

        % ---------------- 分簇布点 ----------------
        case "cluster"
            % 簇数量
            Kc = getfield_def(cfg,'Kc', 6);

            % 簇内扩散尺度（标准差）
            sigma = getfield_def(cfg,'cluster_sigma', min(Lx,Ly)*0.06);

            % 簇中心生成方式
            center_mode = lower(string( ...
                getfield_def(cfg,'cluster_center_mode','random')));

            % 混入均匀节点比例（避免过于理想化）
            mix_u = clamp01( ...
                getfield_def(cfg,'cluster_mix_uniform_ratio', 0.2));

            % 节点数量分配
            Nu = round(N * mix_u);   % 均匀节点
            Nc = N - Nu;             % 簇内节点

            % 生成簇中心
            centers = make_cluster_centers(Kc, Lx, Ly, center_mode);

            % 为每个簇随机分配节点数量
            weights = rand(Kc,1);
            weights = weights / sum(weights);
            sizes = floor(weights * Nc);

            % 修正取整误差
            while sum(sizes) < Nc
                sizes(randi(Kc)) = sizes(randi(Kc)) + 1;
            end
            while sum(sizes) > Nc
                k = find(sizes>0,1);
                sizes(k) = sizes(k)-1;
            end

            % 在每个簇内生成节点
            idx = 1;
            for k = 1:Kc
                nk = sizes(k);
                if nk <= 0, continue; end
                pts = centers(k,:) + sigma * randn(nk,2);
                pos(idx:idx+nk-1,:) = pts;
                idx = idx + nk;
            end

            % 混入均匀分布节点
            if Nu > 0
                pos(idx:end,:) = uniform_points(Nu, Lx, Ly);
            end

            pos = clip_points(pos, Lx, Ly);

            meta.mode = 'cluster';
            meta.cluster_centers = centers;
            meta.cluster_sigma = sigma;

        % ---------------- 沿道路布点 ----------------
        case "road"
            % 均匀节点比例
            mix_u = clamp01( ...
                getfield_def(cfg,'road_mix_uniform_ratio', 0.2));

            Nu = round(N * mix_u);   % 均匀节点
            Nr = N - Nu;             % 非均匀节点

            % 道路宽度（横向高斯扩散）
            road_width = getfield_def(cfg,'road_width', min(Lx,Ly)*0.01);

            % 沿道路分布比例
            road_node_ratio = clamp01( ...
                getfield_def(cfg,'road_node_ratio', 0.8));

            Nroad = round(Nr * road_node_ratio);
            Noff  = Nr - Nroad;

            % 道路定义
            if isfield(cfg,'roads') && ~isempty(cfg.roads)
                roads = cfg.roads;
            else
                road_count = getfield_def(cfg,'road_count', 2);
                roads = generate_random_roads(road_count, Lx, Ly);
            end

            % 沿道路采样
            pts_road = sample_points_along_roads(Nroad, roads, road_width);

            % 非道路节点
            pts_off = uniform_points(Noff, Lx, Ly);

            % 合并
            pos = [pts_road; pts_off; uniform_points(Nu, Lx, Ly)];
            pos = clip_points(pos, Lx, Ly);

            meta.mode = 'road';
            meta.roads = roads;
            meta.road_width = road_width;

        otherwise
            error("未知拓扑模式: %s", topo_mode);
    end

    % =========================================================
    %           热点区域（一个区域很多传感器）
    % =========================================================
    hotspot_enable = getfield_def(cfg,'hotspot_enable', false);
    if hotspot_enable
        ratio = clamp01(getfield_def(cfg,'hotspot_ratio', 0.2));
        Nh = round(N * ratio);
        sigma = getfield_def(cfg,'hotspot_sigma', min(Lx,Ly)*0.05);

        if isfield(cfg,'hotspot_center') && ~isempty(cfg.hotspot_center)
            c = cfg.hotspot_center(:)';
        elseif isfield(cfg,'fire_pos') && ~isempty(cfg.fire_pos)
            c = cfg.fire_pos(:)';
        else
            c = [Lx/2, Ly/2];
        end

        % 将部分节点重新分布到热点区域
        ids = randperm(N, Nh);
        pos(ids,:) = c + sigma * randn(Nh,2);
        pos = clip_points(pos, Lx, Ly);

        meta.hotspot_enable = true;
        meta.hotspot_center = c;
        meta.hotspot_ratio = ratio;
    end

    % =========================================================
    %                  输出结果
    % =========================================================
    topo.pos_sensors = pos;
    topo.BS_pos = BS_pos;
    topo.pos_all = [pos; BS_pos];   % 最后一行为基站
    topo.meta = meta;
end
% ===== 本文件本地函数：并行 worker 可用 =====
function mustHave(s, names)
    for i=1:numel(names)
        if ~isfield(s, names{i})
            error("generate_topology: missing required field '%s'", names{i});
        end
    end
end

function v = getfield_def(s, name, def)
    if isfield(s,name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = def;
    end
end

function x = clamp01(x)
    x = max(0, min(1, x));
end
function pos = uniform_points(N, Lx, Ly)
    pos = [Lx*rand(N,1), Ly*rand(N,1)];
end

function pos = clip_points(pos, Lx, Ly)
    pos(:,1) = max(0, min(Lx, pos(:,1)));
    pos(:,2) = max(0, min(Ly, pos(:,2)));
end

function centers = make_cluster_centers(K, Lx, Ly, mode)
    switch mode
        case "random"
            centers = [Lx*rand(K,1), Ly*rand(K,1)];
        case "grid"
            g = ceil(sqrt(K));
            xs = linspace(0.15*Lx, 0.85*Lx, g);
            ys = linspace(0.15*Ly, 0.85*Ly, g);
            [X,Y] = meshgrid(xs,ys);
            C = [X(:), Y(:)];
            centers = C(1:K,:);
        otherwise
            error("Unknown cluster_center_mode: %s", mode);
    end
end

function roads = generate_random_roads(road_count, Lx, Ly)
    roads = cell(road_count,1);
    for r=1:road_count
        if rand < 0.5
            y0 = Ly*(0.2 + 0.6*rand);
            x1 = Lx*(0.05 + 0.1*rand);
            x2 = Lx*(0.85 + 0.1*rand);
            y1 = y0 + Ly*0.1*(randn);
            y2 = y0 + Ly*0.1*(randn);
            roads{r} = [x1 y1; (x1+x2)/2 y0; x2 y2];
        else
            x0 = Lx*(0.2 + 0.6*rand);
            y1 = Ly*(0.05 + 0.1*rand);
            y2 = Ly*(0.85 + 0.1*rand);
            x1 = x0 + Lx*0.1*(randn);
            x2 = x0 + Lx*0.1*(randn);
            roads{r} = [x1 y1; x0 (y1+y2)/2; x2 y2];
        end
        roads{r}(:,1) = max(0,min(Lx,roads{r}(:,1)));
        roads{r}(:,2) = max(0,min(Ly,roads{r}(:,2)));
    end
end

function pts = sample_points_along_roads(N, roads, road_width)
    if N <= 0
        pts = zeros(0,2);
        return;
    end

    segs = []; % [road_id, p1x,p1y,p2x,p2y, len, nx, ny]
    for rid=1:numel(roads)
        poly = roads{rid};
        for k=1:size(poly,1)-1
            p1 = poly(k,:);
            p2 = poly(k+1,:);
            d = p2 - p1;
            len = hypot(d(1), d(2));
            if len < 1e-9, continue; end
            n = [-d(2), d(1)]/len; % 法向量
            segs = [segs; rid, p1, p2, len, n]; %#ok<AGROW>
        end
    end

    if isempty(segs)
        pts = zeros(N,2);
        return;
    end

    lens = segs(:,6);
    cdf = cumsum(lens) / sum(lens);

    pts = zeros(N,2);
    for i=1:N
        u = rand;
        sid = find(cdf >= u, 1, 'first');
        p1 = segs(sid,2:3);
        p2 = segs(sid,4:5);
        n  = segs(sid,7:8);

        t = rand;
        p = p1 + t*(p2 - p1);

        off = road_width * randn;
        pts(i,:) = p + off*n;
    end
end