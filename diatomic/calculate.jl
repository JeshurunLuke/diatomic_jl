using LinearAlgebra
using Tullio
using Base.Threads
function transition_dipole_moment(H::MoleculeHamiltonian, gs, eigenstates; M = [-1, 0, 1], comp = [1, 1, 1])
    dipoleMat = DipoleMatrix(H.MolOp)
    tdm = zeros(Float64, length(gs), 3)
    println(size(dipoleMat[1]), " ", size(eigenstates),  " ", size(dipoleMat[1]*eigenstates))
    @threads for m_i in M
        d = dipoleMat[m_i + 2] 
        @tullio tdm_t[k] := gs[i]*d[i, j]*eigenstates[j, k] #eigenstates*dipoleMat[m_i + 2])*gs#[dot(gs, temp[:, i]) for i in 1:size(temp, 2)] #dot(gs, dipoleMat[m_i + 2]*eigenstates)
        tdm[:, m_i + 2] = real.(tdm_t)*comp[m_i + 2]
    end
    tdm
end


function transition_dipole_moment(H::MoleculeHamiltonian, gs, eigenstates,  M)
    dipoleMat = DipoleMatrix(H.MolOp)
    tdm = zeros(Float64, length(gs), 3)
    d = dipoleMat[M + 2]
    @tullio tdm_t[k] := gs[i]*d[i, j]*eigenstates[j, k]
end


function magnetic_moment(H::MoleculeHamiltonian, eigenstates)
    muz = zeeman_ham(H.MolOp, [0, 0, 1.0])
    @tullio mu[k] := conj.(eigenstates)[j, k]*muz[j, l]*eigenstates[l, k]
end

function electric_moment(H::MoleculeHamiltonian, eigenstates)
    dz = -1*dc(H.MolOp, [0, 0, 1.0])
    @tullio d[k] := conj.(eigenstates)[j, k]*dz[j, l]*eigenstates[l, k]
end




function findState(stateOI::State, eigsol::sol)
    basisState = Ket(stateOI)
    findMaxOverlap(basisState, eigsol.vec)
end

function findState(width::Vector{T}, eigenergy::Vector{T}) where {T<:Real}
    findall(eigenergy .> minimum(width) .&& eigenergy .< maximum(width))
end



function findTransition(H::MoleculeHamiltonian, stateOI::State, eigsol::sol, TransitionEnergy; M = [-1, 0, 1], width = 5e3, strength = 10, comp = [1, 1, 1])
    indOI = findState(stateOI, eigsol)

    println("Starting from $(KetName(eigsol.vec[:, indOI], basisUC))")
    tdm = transition_dipole_moment(H, eigsol.vec[:, indOI], eigsol.vec, M = M, comp = normalize!(comp))


    energyCenter = real.(eigsol.val[indOI] + TransitionEnergy)
    indstest = findState([energyCenter - width, energyCenter + width], real.(eigsol.val))
    lim = maximum(sum(abs.(tdm[indstest, :]), dims = 2))
    indstest = indstest[ vec(sum(abs.(tdm[indstest, :]), dims = 2)) .>= lim/strength]
    
    println("Found $(length(indstest))")
    indstest, eigsol.vec[:, indstest], eigsol.val[indstest] .-eigsol.val[indOI], tdm[indstest, :] 
end

function findMaxOverlap(basisState, eigvec)
    overlap = 0
    indOI = 0 
    for ind in 1:size(eigvec, 2)
        overlapcurr = abs2.(dot(basisState, eigvec[:, ind]))
        if overlapcurr >= overlap
            overlap = overlapcurr
            indOI = ind
        end
    end
    indOI
end