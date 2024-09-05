
function IsingHam(L::Int64;J::Number=1,h::Number=0,hz::Number = 1e-2)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    σx = [0 1;1 0]
    σz = [1 0;0 -1]

    d = 2
    D = 3

    MPO = Vector{AbstractTensorMap}(undef, L)

    phys = ℂ^d
    bond = ℂ^D
    for i in 1:L
        if i == 1
            #M = TensorMap(vcat(map(x -> x[:],(h*σx,J*σz,I))...), phys' → phys'  ⊗ bond)
            M = TensorMap(vcat(map(x -> x[:],(h*σx+hz*σz,J*σz,I))...), bond' → phys' ⊗ phys)
            M = permute(M,(1,3),(2,))
        elseif i == L
            #M = TensorMap(vcat(map(x -> x[:],(I,σz,σx))...), phys' ⊗ bond → phys' )
            M = TensorMap(vcat(map(x -> x[:],(I,σz,h*σx+hz*σz))...), bond → phys ⊗ phys' )
            M = permute(M,(2,),(1,3))
        else
            Hi = [
                I I0 I0
                σz I0 I0
                h*σx+hz*σz J*σz I
            ]
            reshape(Hi[:],D,d,D,d)
            M = TensorMap(Hi,phys' ⊗ bond'→ phys' ⊗ bond' )
            M = permute(M,(1,4),(3,2))
        end
        MPO[i] = M
        println("MPO finished $i/$L")
    end
    println("MPO totally finished")

    return MPO
end

function IsingMagmom(L::Int64)
    
    I = diagm(ones(2))
    I0 = zeros(2,2)
    σz = [1 0;0 -1]

    d = 2
    D_MPO = 2

    MPO = Vector{AbstractTensorMap}(undef, L)

    phys = ℂ^d
    bond = ℂ^D_MPO
    for i in 1:L
        if i == 1
            #M = TensorMap(vcat(map(x -> x[:],(σz,I))...), phys' → phys'  ⊗ bond) final format
            M = TensorMap(vcat(map(x -> x[:],(σz,I))...), bond' → phys'  ⊗ phys)
            M = permute(M,(1,3),(2,))
        elseif i == L
            #M = TensorMap(vcat(map(x -> x[:],(I,σz))...), phys' ⊗ bond → phys' ) final format
            M = TensorMap(vcat(map(x -> x[:],(I,σz))...), bond → phys ⊗ phys' )
            M = permute(M,(2,),(1,3))
        else
            Hi = [
                I I0
                σz I
            ]
            reshape(Hi[:],D_MPO,d,D_MPO,d)
            M = TensorMap(Hi,phys' ⊗ bond'→ phys' ⊗ bond' )
            M = permute(M,(1,4),(3,2))
        end
        MPO[i] = M
        println("MPO finished $i/$L")
    end
    println("MPO totally finished")

    return MPO

end


function IsingMPS(L::Int64,state::String,D_MPS::Int64;noise::Number=0)

    initialState = zeros(ComplexF64,[d for _ in 1:L]...)

    if state == "FM"
        initialState[[1 for _ in 1:L]...] = 1im
    elseif state == "AFM"
        initialState[repeat([1, 2], div(L,2))...] = 1im
    else
        @error "state not defined"
    end

    initialState += noise .* rand(ComplexF64,[d for _ in 1:L]...)

    MPS = let 
        iniMPS = Tensor(initialState, ⊗([ℂ^d for i in 1:L]...)) |> x -> x / norm(x)
        MPS = Vector{AbstractTensorMap}(undef, L)
        for ii in L:-1:2
            if ii == L
                U,S,V = tsvd(iniMPS,Tuple.((1:ii-1,ii:L))...;trunc = truncdim(D_MPS))
                MPS[ii] = V
            else
                U,S,V = tsvd(iniMPS,Tuple.((1:ii-1,ii:ii+1))...;trunc = truncdim(D_MPS))
                MPS[ii] = V
            end
            iniMPS = U*S
            println("MPS initialized $(L-ii+1)/$L")
        end
        MPS[1] = permute(iniMPS,(),(1,2))
        println("MPS initialized $L/$L")
        println("MPS totally initialized")
        MPS
    end

    return MPS
end



