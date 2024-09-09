using LinearAlgebra

# 参数定义
t = 1.0  # 跃迁能量
a = 1.0  # 晶格常数
N = 3000  # 更高的波矢点密度

# 生成二维波矢空间
k_vals = range(-pi/a, stop=pi/a, length=N)
grid_kx, grid_ky = [k_vals for _ in 1:N], [k_vals for _ in 1:N]

# 定义能量色散关系
function energy(kx, ky, t, a)
    return 2t * (cos(kx*a) + cos(ky*a))
end

# 计算能量矩阵
E = [energy(kx, ky, t, a) for kx in k_vals, ky in k_vals]

# 统计态密度
function density_of_states(E, target_energy, deltaE=0.001)
    count = sum(1 for e in E if abs(e - target_energy) < deltaE)
    return count / (N^2)  # 归一化
end

# 计算在Gamma点的态密度
E_gamma = energy(0, 0, t, a)
dos_at_gamma = density_of_states(E, E_gamma)

println("态密度在 Gamma 点 (数值计算): $dos_at_gamma")
