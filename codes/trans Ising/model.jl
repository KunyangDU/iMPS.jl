
function tranIsingH(L::Int64,d::Int64,D::Int64;J::Number=1,h::Number=0)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    σx = [0 1;1 0]
    σz = [1 0;0 -1]

    MPO = Vector{AbstractTensorMap}(undef, L)

    phys = ℂ^d
    bond = ℂ^D
    for i in 1:L
        if i == 1
            M = TensorMap(vcat(map(x -> x[:],(h*σx,J*σz,I))...), phys' → phys'  ⊗ bond)
        elseif i == L
            M = TensorMap(vcat(map(x -> x[:],(I,σz,σx))...), phys' ⊗ bond → phys' )
        else
            Hi = [
                I I0 I0
                σz I0 I0
                h*σx J*σz I
            ]
            reshape(Hi[:],D,d,D,d)
            M = TensorMap(Hi,phys' ⊗ bond→ phys' ⊗ bond )
        end
        MPO[i] = M
        println("MPO finished $i/$L")
    end
    println("MPO totally finished")

    return MPO
end

function IsingMagmom(L::Int64,d::Int64,D_MPO::Int64)
    
    I = diagm(ones(2))
    I0 = zeros(2,2)
    σz = [1 0;0 -1]

    MPO = Vector{AbstractTensorMap}(undef, L)

    phys = ℂ^d
    bond = ℂ^D_MPO
    for i in 1:L
        if i == 1
            M = TensorMap(vcat(map(x -> x[:],(σz,I))...), phys' → phys'  ⊗ bond)
        elseif i == L
            M = TensorMap(vcat(map(x -> x[:],(I,σz))...), phys' ⊗ bond → phys' )
        else
            Hi = [
                I I0
                σz I
            ]
            reshape(Hi[:],D_MPO,d,D_MPO,d)
            M = TensorMap(Hi,phys' ⊗ bond→ phys' ⊗ bond )
        end
        MPO[i] = M
        println("MPO finished $i/$L")
    end
    println("MPO totally finished")

    return MPO

end



