import Base.+
import Base.*
mutable struct State 
    comp::Array{ComplexF64, 2}
    QN::Array{Float64, 2}
    basis::Array{Float64, 2}
end
function *(b::T, a::State) where T <: Number
    State(b*a.comp, a.QN, a.basis)
end
function norm(a::State)
    sum(abs2.(a.comp))
end
function Ket(state::State) 
    Kets = zeros(ComplexF64, size(state.basis, 1))
    indsOI = [findall(row -> all(row .== state.QN[comp_ind, :]), eachrow(state.basis)) for comp_ind in 1:size(state.QN, 1)]
    indsOI = vcat(indsOI...)
    Kets[indsOI] = state.comp
    return vec(Kets)
end
function +(a::State ,  b::State) 
    if !(length(a.basis) == length(b.basis))
        return "Error: Mismatch in length"
    end
    

    indsOI_a = indexin([vec(i) for i in eachrow(a.QN)], [vec(i) for i in eachrow(a.basis)])
    comp_a = zeros(ComplexF64, size(a.basis)[1])
    comp_a[indsOI_a] = a.comp
    indsOI_b = indexin([vec(i) for i in eachrow(b.QN)], [vec(i) for i in eachrow(a.basis)])
    comp_b = zeros(ComplexF64, size(a.basis)[1])
    comp_b[indsOI_b] = b.comp

    return State(hcat((comp_a .+ comp_b)...) , a.basis, a.basis)
end
function KetName(state::Vector{<:Any}, basis; QMorder = [5, 3, 2, 1])
    indsorted = sortperm(abs.(state), rev = true)
    
    prod_str = join([string("$(round(abs2.(state[indsorted[ind]]), digits = 3)) | ", join(basis[indsorted[ind]][QMorder], ", "), " >") for ind in 1:3], " + ")
    
end

State(comp::Vector{<:Any}, QN::Vector{<:Any},     basis::Vector{<:Any}) = State(ComplexF64.(hcat(comp...)), Matrix(hcat(QN...)'), Matrix(hcat(basis...)'))


