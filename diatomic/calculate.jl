using LinearAlgebra
using Tullio
using Base.Threads
function transition_dipole_moment(H::MoleculeHamiltonian, gs, eigenstates; M = [-1, 0, 1], comp = [1, 1, 1])
    dipoleMat = DipoleMatrix(H.MolOp)
    tdm = zeros(Float64, length(gs), 3)
    
    @threads for m_i in M
        d = dipoleMat[m_i + 2] 
        @tullio tdm_t[k] := gs[i]*d[i, j]*eigenstates[j, k] #eigenstates*dipoleMat[m_i + 2])*gs#[dot(gs, temp[:, i]) for i in 1:size(temp, 2)] #dot(gs, dipoleMat[m_i + 2]*eigenstates)
        tdm[:, m_i + 2] = real.(tdm_t)*comp[m_i + 2]
    end
    tdm
end



label_states_I_MI(Hmol::MoleculeHamiltonian, state::State) = State(label_states_I_MI(Hmol, Ket(state)), getBasisFC(Hmol.MolOp.basisTree,"I_c"), getBasisFC(Hmol.MolOp.basisTree,"I_c"))
label_states_I_MI(Hmol::MoleculeHamiltonian, eigsol::sol) = label_states_I_MI(Hmol, eigsol.vec)
function label_states_I_MI(Hmol::MoleculeHamiltonian, state::Vector{<:ComplexF64})
    UnitaryUC2I = step_up(Hmol, "I_c")
    UnitaryUC2I*state
end

function label_states_I_MI(Hmol::MoleculeHamiltonian, eigstate::Matrix{<:ComplexF64})
    UnitaryUC2I = step_up(Hmol, "I_c")
    UnitaryUC2I*eigstate*adjoint(UnitaryUC2I)
end



label_states_F_MF(Hmol::MoleculeHamiltonian, state::State) = State(label_states_F_MF(Hmol, Ket(state)), getBasisFC(Hmol.MolOp.basisTree,"F"), getBasisFC(Hmol.MolOp.basisTree,"F"))
label_states_F_MF(Hmol::MoleculeHamiltonian, eigsol::sol) = label_states_F_MF(Hmol, eigsol.vec)
function label_states_F_MF(Hmol::MoleculeHamiltonian, state::Vector{<:ComplexF64})
    UnitaryUC2F = step_up(Hmol, "F")
    UnitaryUC2F*state
end

function label_states_F_MF(Hmol::MoleculeHamiltonian, eigstate::Matrix{<:ComplexF64})
    UnitaryUC2F = step_up(Hmol, "F")
    UnitaryUC2F*eigstate*adjoint(UnitaryUC2F)
end



function transition_dipole_moment(H::MoleculeHamiltonian, gs, eigenstates,  M)
    dipoleMat = DipoleMatrix(H.MolOp)
    tdm = zeros(Float64, length(gs), 3)
    d = dipoleMat[M + 2]
    @tullio tdm_t[k] := gs[i]*d[i, j]*eigenstates[j, k]
end





function magnetic_moment(H::MoleculeHamiltonian, eigenstates)
    muz = H.Hzeem#zeeman_ham(H.MolOp, [0, 0, 1.0])
    @tullio mu[k] := conj.(eigenstates)[j, k]*muz[j, l]*eigenstates[l, k]
end

function electric_moment(H::MoleculeHamiltonian, eigenstates)
    dz = H.Hdc#-1*dc(H.MolOp, [0, 0, 1.0])
    @tullio d[k] := conj.(eigenstates)[j, k]*dz[j, l]*eigenstates[l, k]
end


function findState(stateOI::Vector{<:Complex}, eigsol::sol)
    findMaxOverlap(stateOI, eigsol.vec)
end

function findState(stateOI::State, eigsol::sol)
    basisState = Ket(stateOI)
    findMaxOverlap(basisState, eigsol.vec)
end
function findState(width::Vector{T}, eigenergy::Vector{T}) where {T<:Real}
    findall(eigenergy .> minimum(width) .&& eigenergy .< maximum(width))
end

function findTransition(H::MoleculeHamiltonian, stateOI::State, eigsol::sol, TransitionEnergy::Number; M = [-1, 0, 1], width = 5e3, strength = 10, comp = [1, 1, 1])
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

findTransition(H::MoleculeHamiltonian, stateOI::State, eigsol::sol;  M = [-1, 0, 1], NumOfStates = 10, comp = [1.0, 1, 1]) = findTransition(H, Ket(stateOI), eigsol;  M = M, NumOfStates = NumOfStates, comp = comp)
findTransition(H::MoleculeHamiltonian, stateOI::State, eigsol::sol, N::Vector{<:Int};  M = [-1, 0, 1], NumOfStates = 10, comp = [1.0, 1, 1]) = findTransition(H, Ket(stateOI), eigsol, N;  M = M, NumOfStates = NumOfStates, comp = comp)
function findTransition(H::MoleculeHamiltonian, stateOI::Vector{<:Complex}, eigsol::sol;  M = [-1, 0, 1], NumOfStates = 50, comp = [1.0, 1, 1], display = true)
    indOI = findState(stateOI, eigsol)
    basisUC = getBasisUC(H.MolOp.basisTree)
    println("Starting from $(KetName(eigsol.vec[:, indOI], basisUC))")
    tdm = transition_dipole_moment(H, eigsol.vec[:, indOI], eigsol.vec, M = M, comp = normalize!(comp))

    indOrder = reverse(sortperm(vec(sum(abs.(tdm[:, :]), dims = 2))))
    stateName, eigvec, eigval, tdm = indOrder, eigsol.vec[:, indOrder[1:NumOfStates]], eigsol.val[indOrder[1:NumOfStates]] .-eigsol.val[indOI], tdm[indOrder[1:NumOfStates], :] 


    if display
        for state in 1:length(eigval)
            println("State: ", KetName(eigvec[:, state], basisUC), "| Energy (MHz): ", eigval[state]*1e-6)
            println("\t TDM [sigma_m, pi, sigma_p] = $(tdm[state, :])")
        end
    end


    stateName, eigvec, eigval, tdm
end




function findTransition(H::MoleculeHamiltonian, stateOI::Vector{<:Complex}, eigsol::sol, N::Vector{<:Int};  M = [-1, 0, 1], NumOfStates = 50, comp = [1.0, 1, 1], display = true)
    indOI = findState(stateOI, eigsol)
    basisUC = getBasisUC(H.MolOp.basisTree)
    println("Starting from $(KetName(eigsol.vec[:, indOI], basisUC))")
    tdm = transition_dipole_moment(H, eigsol.vec[:, indOI], eigsol.vec, M = M, comp = normalize!(comp))

    spinDim = (prod([length(NuclearSpin.spin) for NuclearSpin in endNode(H.MolOp.basisTree)[2:end]]))
    NOI = vcat([[(sum(spinDim.*[(2*N_i2 + 1) for N_i2 in 0:(N_i-1)]) + 1):sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:N_i])...] for N_i in N]...)
    
    indOrder = reverse(sortperm(vec(sum(abs.(tdm[:, :]), dims = 2))))
    indOrder = intersect(indOrder, NOI)

    stateName, eigvec, eigval, tdm = indOrder, eigsol.vec[:, indOrder[1:NumOfStates]], eigsol.val[indOrder[1:NumOfStates]] .-eigsol.val[indOI], tdm[indOrder[1:NumOfStates], :] 

    if display
        for state in 1:length(eigval)
            println("State: ", KetName(eigvec[:, state], basisUC), "\t| Energy (MHz): ", eigval[state]*1e-6)
            println("\t TDM [sigma_m, pi, sigma_p] = $(tdm[state, :])")
        end
    end

    stateName, eigvec, eigval, tdm
end






diabaticRamp(startingState::State, eigsol_vec::Vector{sol}, Field_ramp::Vector{Float64}; Field = "B") = diabaticRamp(Ket(startingState), eigsol_vec, Field_ramp, Field = Field)
diabaticRamp(Hmol::MoleculeHamiltonian, startingState::State, eigsol_vec::Vector{sol}, Field_ramp::Vector{Float64}; Field = "B") = State(diabaticRamp(Ket(startingState), eigsol_vec, Field_ramp, Field = Field), Hmol)
function diabaticRamp(startingState::Vector{<:Complex}, eigsol_vec::Vector{sol}, Field_ramp::Vector{Float64}; Field = "B")
    Field_s = zeros(Float64, length([1 for i in eigsol_vec]))#Vector{Float64}[]
    if Field == "B"
        Field_s = [eigsol.B_field for eigsol in eigsol_vec]
    else Field == "E"
        Field_s = [eigsol.E_field for eigsol in eigsol_vec]
    end

        #[eigsol.B_Field for eigsol in eigsol_vec]
    indStart = argmin(abs.(Field_s .- Field_ramp[1]))
    indEnd = argmin(abs.(Field_s .- Field_ramp[2]))

    moleculeind = findMaxOverlap(startingState, eigsol_vec[indStart].vec)
    for ind_c in (indStart-1):-1:indEnd
        moleculeind = findMaxOverlap(eigsol_vec[ind_c + 1].vec[:, moleculeind], eigsol_vec[ind_c].vec)#getMaxOverlap(stateOI_n, composition_m[ind, :, :])
    end

    eigsol_vec[indEnd].vec[:, moleculeind]
end



scanDiabiaticEnergy(startingState::State, eigsol_vec::Vector{sol}, Field_ramp::Vector{Float64}; Field = "B") = scanDiabiaticEnergy(Ket(startingState), eigsol_vec, Field_ramp, Field = Field)
scanDiabiaticEnergy(Hmol::MoleculeHamiltonian, startingState::State, eigsol_vec::Vector{sol}, Field_ramp::Vector{Float64}; Field = "B") = scanDiabiaticEnergy(Ket(startingState), eigsol_vec, Field_ramp, Field = Field)

function scanDiabiaticEnergy(startingState::Vector{<:Complex}, eigsol_vec::Vector{sol}, Field_ramp::Vector{Float64}; Field = "B")
    Field_s = zeros(Float64, length([1 for i in eigsol_vec]))#Vector{Float64}[]
    if Field == "B"
        Field_s = [eigsol.B_field for eigsol in eigsol_vec]
    else Field == "E"
        Field_s = [eigsol.E_field for eigsol in eigsol_vec]
    end

        #[eigsol.B_Field for eigsol in eigsol_vec]
    indStart = argmin(abs.(Field_s .- Field_ramp[1]))
    indEnd = argmin(abs.(Field_s .- Field_ramp[2]))
    energyList = []
    moleculeind = findMaxOverlap(startingState, eigsol_vec[indStart].vec)
    push!(energyList, eigsol_vec[indStart].val[moleculeind])
    for ind_c in (indStart-1):-1:indEnd
        moleculeind = findMaxOverlap(eigsol_vec[ind_c + 1].vec[:, moleculeind], eigsol_vec[ind_c].vec)#getMaxOverlap(stateOI_n, composition_m[ind, :, :])
        push!(energyList, eigsol_vec[ind_c].val[moleculeind])
    end
    reverse(energyList)
end



findAvoidedCrossing(startingState::State, eigsol_vec::Vector{sol}, Field_ramp::Vector{Float64}; Field = "B") = findAvoidedCrossing(Ket(startingState), eigsol_vec, Field_ramp, Field = Field)

function findAvoidedCrossing(startingState::Vector{<:Complex}, eigsol_vec::Vector{sol}, Field_ramp::Vector{Float64}; Field = "B")
    Field_s = zeros(Float64, length([1 for i in eigsol_vec]))
    if Field == "B"
        Field_s = [eigsol.B_field for eigsol in eigsol_vec]
    else Field == "E"
        Field_s = [eigsol.E_field for eigsol in eigsol_vec]
    end

        #[eigsol.B_Field for eigsol in eigsol_vec]
    indStart = argmin(abs.(Field_s .- Field_ramp[1]))
    indEnd = argmin(abs.(Field_s .- Field_ramp[2]))

    moleculeind = findMaxOverlap(startingState, eigsol_vec[indStart].vec)
    avoided_crossing = []
    for ind_c in (indStart - 1):-1:indEnd
        moleculeind_c = findMaxOverlap(eigsol_vec[ind_c + 1].vec[:, moleculeind], eigsol_vec[ind_c].vec)#getMaxOverlap(stateOI_n, composition_m[ind, :, :])
        if abs.(moleculeind_c - moleculeind) > 0
            energyGap =  abs.(eigsol_vec[ind_c + 1].val[moleculeind] - eigsol_vec[ind_c].val[moleculeind_c])
            push!(avoided_crossing, [(Field_s[ind_c] + Field_s[ind_c + 1])/2,energyGap])
        end
        moleculeind = moleculeind_c
    end
    avoided_crossing
end