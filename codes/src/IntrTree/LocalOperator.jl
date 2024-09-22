abstract type AbstractLocalOperator end

mutable struct IdentityOperator <: AbstractLocalOperator
    Opri::Union{Nothing, AbstractTensorMap}
    site::Int64
    strength::Number 
    function IdentityOperator(Opri::AbstractTensorMap,site::Int64, strength::Number = NaN)
        return new(Opri, site, strength)
   end
    function IdentityOperator(site::Int64, strength::Number = NaN)
         return new(nothing, site, strength)
    end
end

function Base.show(io::IO,Opr::IdentityOperator)
    print(io,"I$(String(collect("$(Opr.site)") .+ 8272))")
    if !isnan(Opr.strength)
        print(io, "($(Opr.strength))")
   end
end

mutable struct LocalOperator <: AbstractLocalOperator
    Opri::Union{Nothing,AbstractTensorMap}
    name::String
    site::Int64
    strength::Number

    function LocalOperator(Opri::AbstractTensorMap, name::String, site::Int64)
        return new(Opri, name, site, NaN)
    end

    function LocalOperator(Opri::AbstractTensorMap, name::String, site::Int64, strength::Number)
        return new(Opri, name, site, strength)
    end
end

function Base.show(io::IO,Opr::LocalOperator)
    print(io,"$(Opr.name)$(String(collect("$(Opr.site)") .+ 8272))")
    if !isnan(Opr.strength)
        print(io, "($(Opr.strength))")
   end
end


isequal(::AbstractLocalOperator, ::AbstractLocalOperator) = false
isequal(A::IdentityOperator, B::IdentityOperator) = (A.site == B.site)
function isequal(A::LocalOperator, B::LocalOperator)
    A.name ≠ B.name && return false
    A.site ≠ B.site && return false
    A.strength ≠ B.site && return false
    return A.Opri == B.Opri
end

function getIdTensor(Opr::LocalOperator)
    space = domain(Opr.Opri)[1]
    return TensorMap(diagm(ones(dim(space))),space,space)
end
