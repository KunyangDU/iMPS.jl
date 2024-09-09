using LinearAlgebra

# Krylov子空间方法
function Krylov(A, b, ω; max_iter=50, tol=1e-8)
    n = length(b)
    I = diagm(ones(n))
    M = A - ω * I

    # 初始化
    Q = zeros(n, max_iter + 1)
    H = zeros(max_iter + 1, max_iter)
    Q[:, 1] = b / norm(b)  # 初始向量

    # Arnoldi迭代
    for k in 1:max_iter
        # 计算 M * Q[:, k]
        y = M * Q[:, k]
        
        # 正交化
        for j in 1:k
            H[j, k] = dot(Q[:, j], y)
            y -= H[j, k] * Q[:, j]
        end
        
        H[k + 1, k] = norm(y)
        if H[k + 1, k] < tol
            break
        end
        
        Q[:, k + 1] = y / H[k + 1, k]
    end
    
    # 求解上三角系统 H * y = ||b|| * e1
    e1 = zeros(max_iter + 1)
    e1[1] = norm(b)
    @show e1
    y = H[1:max_iter, 1:max_iter] ./ e1[1:max_iter]

    # 计算最终解 x
    x = Q[:, 1:max_iter] * y
    return x
end

# 示例用法
function main()
    n = 100  # 矩阵的维度
    A = rand(n, n) + diagm(ones(n))  # 随机矩阵加上单位矩阵确保A是非奇异的
    ω = 0.9
    b = ones(n)
    
    # 使用Krylov方法计算 (A - ω I)^-1 * b
    x = Krylov(A, b, ω)

    println("Solution vector x:")
    println(x)
end

main()