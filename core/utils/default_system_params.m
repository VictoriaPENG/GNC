function def = default_system_params()
%DEFAULT_SYSTEM_PARAMS 系统级默认参数（仿真环境 + 信道/能耗 + 业务/队列 + 随机种子等）
% 输出:
%   def : struct，包含仿真所需的所有默认参数

def = struct();

%% ========== 1) 场景与拓扑参数 ==========
def.Lx = 1000; 
def.Ly = 1000;
% Lx, Ly : 仿真区域尺寸（通常单位 m）。节点随机分布/生成拓扑的边界范围。

def.N  = 700;
% N : 传感器节点数量（不含基站）。影响网络密度、连通性、拥塞与能耗。

def.Rc = 70;
% Rc : 通信半径/邻居判定半径（m）。用于构建邻接关系、路由可达性等。

def.T  = 2000;
% T : 仿真时长（以"离散时隙/步数"为单位）。每一步通常代表一次业务产生/转发周期。

%% ========== 2) 基站与火源/热点相关参数 ==========
def.BS_pos   = [500, 500];
% BS_pos : 基站坐标 [x, y]（m）。所有数据最终要汇聚到该位置。

def.fire_pos = [250, 250];
% fire_pos : 火源/事件中心坐标 [x, y]（m）。用于建模事件区域或风险区域。

def.Rf = 220;
% Rf : 火源/事件影响半径（m）。半径内节点可能产生更高业务、或受到更高损毁概率等。

%% ========== 2.1) 拓扑生成与热点/分簇/道路参数 ==========
def.topo_mode = 'uniform';
% topo_mode : 拓扑生成模式
%   - 'uniform'（随机均匀）/ 'cluster'（分簇）/ 'road'（沿道路）

% hotspot：支持"一个区域很多传感器"（可叠加到任意 topo_mode）
def.hotspot_enable = false;
def.hotspot_ratio  = 0.25;                   % 热点节点比例
def.hotspot_sigma  = min(def.Lx,def.Ly)*0.05;% 热点扩散尺度（标准差）
def.hotspot_center = def.fire_pos;           % 默认以火区为热点中心

% cluster 参数
def.Kc = 6;                                  % 簇数量
def.cluster_sigma = min(def.Lx,def.Ly)*0.06; % 簇内扩散尺度
def.cluster_center_mode = 'random';          % 'random' 或 'grid'
def.cluster_mix_uniform_ratio = 0.2;         % 混入均匀节点比例

% road 参数
def.road_count = 2;                          % 自动生成道路数量（未提供 roads 时）
def.road_width = min(def.Lx,def.Ly)*0.01;    % 道路横向扩散尺度（标准差）
def.road_mix_uniform_ratio = 0.2;            % 混入均匀节点比例
def.road_node_ratio = 0.8;                   % 非均匀部分中，沿道路采样比例

%% ========== 3) 业务产生参数 ==========
def.lambda = 0.05;
% lambda : 节点每时隙产生数据包的概率/到达率（常见建模为 Bernoulli 或 Poisson 强度）。
%          值越大 -> 业务越繁忙 -> 队列更容易拥塞，时延/丢包上升。

%% ========== 4) 链路质量/PRR 模型参数 ==========
def.alpha = 0.01;
% alpha : PRR 随距离衰减的系数（通常用于 exp(-alpha*d) 或类似模型）。
%         alpha 越大 -> 距离稍远 PRR 就下降更快。

def.prr_floor = 0.08;
% prr_floor : PRR 下限（地板值）。避免出现 0 导致 ETX 无穷大或网络完全不可用。

def.noise_sigma = 0.02;
% noise_sigma : PRR 噪声标准差（随机遮挡/干扰等不确定性）。用于模拟链路抖动。

%% ========== 5) 邻接/可用链路阈值 ==========
def.USE_PRR_TH = true;
% USE_PRR_TH : 是否启用 PRR 阈值过滤邻居。
%              true  -> 只保留 PRR 足够好的链路，拓扑更"干净"但可能更稀疏；
%              false -> 仅按距离 Rc 建邻接（或由其他规则决定）。

def.PRR_MIN = 0.15;
% PRR_MIN : 最小可接受 PRR 阈值（当 USE_PRR_TH=true 时生效）。
%           小于该值的链路视为不可用/不建边，减少低质量转发带来的重传与能耗。

%% ========== 6) 能量模型参数 ==========
def.E0    = 100;
% E0 : 节点初始能量（归一化单位或 Joule，取决于你系统的能耗标定）。

def.Etx   = 0.9;
% Etx : 每次发送一个包消耗的能量（含一次发送动作，不含可能的重传额外消耗）。

def.Erx   = 0.4;
% Erx : 每次接收一个包消耗的能量（转发场景里，中继节点既接收又发送）。

def.Eidle = 0.001;
% Eidle : 每时隙空闲/待机消耗的能量（维持工作状态的背景耗电）。

%% ========== 7) 节点随机故障/火灾致死模型参数 ==========
def.p_rand = 1e-5;
% p_rand : 节点随机故障概率（每时隙）。用于模拟硬件故障、电量异常、环境因素等。

def.beta_fire = 0.05;
% beta_fire : 火灾影响强度/随时间增长系数（具体含义取决于你后续 fire 模型的实现方式）。
%             常用于控制"火灾造成损坏概率"的增长速度。

def.fire_kill_scale = 1e-4;
% fire_kill_scale : 火灾致死概率的尺度系数（把距离/时间/强度映射到实际击杀概率的倍率）。
%                   值越大 -> 火灾区域内节点更快失效。

%% ========== 8) 队列/转发与 MAC 重传参数 ==========
def.Qmax   = 50;
% Qmax : 每个节点队列最大长度（FIFO）。超过则丢包（拥塞丢弃）。

def.TTLmax = 25;
% TTLmax : 包的最大跳数/生存期（Time-To-Live）。避免环路或过长路径导致无限转发。

def.Ksend = 20;
% Ksend : 每时隙每个节点允许发送的最大包数（服务能力/带宽限制的简化）。
%         值越小 -> 更容易排队；值越大 -> 更"理想带宽"。

def.MAX_RETX = 5;
% MAX_RETX : MAC 层最大重传次数（链路失败后的重试上限）。
%            重传越多 -> PDR 可能上升，但能耗与时延也会增加。

def.RETX_TO_TAIL = true;
% RETX_TO_TAIL : 重传包是否放队列尾部（true 更避免"卡死"）。

%% ========== 8.1) 路由/建树与能量感知参数 ==========
def.route_mode = 'mix';
% route_mode : 路由策略（如 mix/tree/tree_minhop/tree_energy/tree_energy_aodv 等）

def.tree_cost = 'etx';
% tree_cost : 汇聚树代价（'etx' | 'hop'）

def.tree_rebuild_period = inf;
% tree_rebuild_period : 汇聚树重建周期（inf=不重建）

def.tree_energy_aware = false;
% tree_energy_aware : 是否启用能量感知建树（Energy-aware Tree）

def.tree_energy_beta  = 2.0;
% tree_energy_beta : 能量惩罚强度（越大越偏向选择高能量父节点）

def.tree_energy_min_frac = 0.10;
% tree_energy_min_frac : 残余能量归一化下限（避免权重过大）

def.local_switch_enable = false;
% local_switch_enable : 是否启用局部切换/修复（AODV hybrid 模式用）

def.local_switch_prr_min = 0.20;
% local_switch_prr_min : 局部切换最小 PRR 阈值

def.local_switch_energy_min = 0.15;
% local_switch_energy_min : 局部切换最小能量阈值（剩余能量比例）

def.local_switch_gain = 0.10;
% local_switch_gain : 局部切换收益阈值

def.local_switch_hop_weight = 0.10;
% local_switch_hop_weight : hop 惩罚权重

def.local_switch_hop_strict = true;
% local_switch_hop_strict : 是否严格限制 hop 变差

%% ========== 8.2) 统计/能耗跟踪 ==========
def.energy_trace = false;
% energy_trace : 是否记录每时隙能耗时间序列（用于绘图/分析）

%% ========== 9) 可复现实验与日志 ==========
def.seed = 1;
% seed : 随机种子。用于可复现拓扑、链路噪声、业务产生、随机故障等随机过程。

def.verbose = false;
% verbose : 是否输出详细日志/调试信息。true 会打印更多过程信息但仿真更慢。

%% ========== 10) 可视化/绘图参数 ==========
def.topo_modes = {};
% topo_modes : 指定需要展示的拓扑模式列表（空=单一模式）

def.show_links = true;
% show_links : 是否显示连边

def.node_size = max(4, min(12, round(6000 / def.N)));
% node_size : 节点绘图大小（随 N 自适应）

def.node_alpha = 0.55;
% node_alpha : 节点透明度

def.node_color = [0.2 0.2 0.2];
def.fire_node_color = [0.85 0.2 0.2];
def.bs_color = [0.1 0.4 0.9];
def.road_color = [0.2 0.6 0.2];
def.hotspot_color = [0.7 0.2 0.7];

def.fire_circle_style = '--r';
def.fire_circle_width = 1.5;

def.node_size_abc = 10;
% node_size_abc : 三联图专用 marker 大小

def.save_path = '';
% save_path : 输出路径（空=不保存）
end
