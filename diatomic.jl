module diatomic_jl


include("QuantumToolbox/QuantumWrapper.jl")
using .QuantumWrapper
import .QuantumWrapper: AM_Toolbox
import .QuantumWrapper.State_Toolbox.State
import .QuantumWrapper.State_Toolbox.KetName

export Molecule, Hamiltonian, solve, plotting, calculate, AM_Toolbox, getBasisUC, State, KetName

function getBasisUC end



module MoleculeTypes
include("diatomic/MolecularConsts.jl")
export Rb87Cs133, K41Cs133, K40Rb87 ,Na23K40, Na23Rb87, K40Rb87T
end



### Initializes the molecule and the Appropriate spin operators
module Molecule
import ..QuantumWrapper.AM_Toolbox: Basis, Node, angularMomentum, couple, endNode, AngularNode
using ..QuantumWrapper.SpinOperator
using LinearAlgebra
using SparseArrays
export moleculeProperties, generatorMolecule
struct moleculeProperties{T<:Node} 
    Constants
    basisTree::T
    N
    I1
    I2
end

function generateMolecule(MoleculeInfo, Nmax)
    I1 = angularMomentum("I1", [MoleculeInfo["I1"]])
    I2 = angularMomentum("I2", [ MoleculeInfo["I2"]])
    N = angularMomentum("N", [0:Nmax...])
    I_c = couple("I_c", I2, I1)
    F = couple("F", N, I_c)
    Nop, I2op, I1op = generateSpinOperator(F)
    moleculeProperties(MoleculeInfo,F, Nop, I1op, I2op)
end 
function generateSpinOperator(F::Node)
    EndNodes = endNode(F)
    Operators = []
    for (i, node) in enumerate(EndNodes)
        dimsA = i == length(EndNodes) ? 1 : prod([length(node.spin) for node in EndNodes[(i+ 1):end]])
        dimsB = i == 1 ? 1 : prod([length(node.spin) for node in EndNodes[(i-1):-1:1]])



        I_a = sparse(1.0I, dimsA, dimsA)
        I_b = sparse(1.0I, dimsB, dimsB)

        push!(Operators, [kron(I_b, kron(x_operator(node), I_a)), kron(I_b, kron(y_operator(node), I_a)), kron(I_b, kron(z_operator(node), I_a)), kron(I_b, kron(raising_operator(node), I_a))])
    end
    Operators
end
end
##### Used to Generate the Hamiltonian with the spin operators ################################
module Hamiltonian
import ..QuantumWrapper.AM_Toolbox: endNode, getBasis#, getBasisUC
import ..QuantumWrapper: AM_Toolbox
export generateHamiltonian, OpticalBeam, MoleculeHamiltonian, zeeman_ham, dc, DipoleMatrix
import ..Molecule: moleculeProperties
include("diatomic/Hamiltonian.jl")
using ..diatomic_jl
import ..diatomic_jl: getBasisUC
getBasisUC(mol::MoleculeHamiltonian) = AM_Toolbox.getBasisUC(mol.MolOp.basisTree)
end

##### Used to solve the Hamiltonian ################################
module solve
export scanZeeman, scanStark, diagonalize, sol, scanIntensity
import ..Hamiltonian: MoleculeHamiltonian
using LinearAlgebra
using Arpack
using SparseArrays


struct sol
    B_field
    E_field
    Intensity
    vec
    val
end

function diagonalize(H::MoleculeHamiltonian{T}, B_field, E_field; Intensity = 1,  nev = 36) where {T<:AbstractMatrix}
    J2Hz = 1.509190311676e33 
    
    vals, vecs = eigen(Hermitian(Matrix(H.Hhfs .+ H.Hzeem*B_field .+ E_field*H.Hdc + sum(Intensity .* H.Hac))))
    sorted_indices = sortperm(vals)
    vals = vals[sorted_indices]
    vecs = vecs[:, sorted_indices]
    sol(B_field, E_field, [Intensity...], vecs, vals*J2Hz)
end
function diagonalize(H::MoleculeHamiltonian{T}, B_field, E_field; Intensity = 1,  nev = 36) where T<:SparseMatrixCSC
    J2Hz = 1.509190311676e33 
    vals, vecs = eigs(H.Hhfs .+ H.Hzeem*B_field .+ E_field*H.Hdc + sum(Intensity .* H.Hac), nev = nev, which=:SM)
    sorted_indices = sortperm(real.(vals))
    vals = vals[sorted_indices]
    vecs = vecs[:, sorted_indices]
    sol(B_field, E_field, [Intensity...], vecs, real.(vals)*J2Hz)
end

function scanIntensity(H::MoleculeHamiltonian, B_Field::T, E_Field::T, Intensity) where T<:Number
    J2Hz = 1.509190311676e33 
    solutions = Vector{sol}(undef, size(Intensity, 1))
    for (i, Int_i) in enumerate(Intensity)
        solutions[i] = diagonalize(H, B_Field, E_Field, Intensity = Int_i)
    end
    solutions
end


function scanZeeman(H::MoleculeHamiltonian, B_field::Vector{<:Real}, E_field::T; Intensity = 1, nev = 36) where T<:Number
    J2Hz = 1.509190311676e33 
    solutions = Vector{sol}(undef, length(B_field))
    for (i, B_i) in enumerate(B_field)
        solutions[i] = diagonalize(H, B_i , E_field, Intensity = 1, nev = nev)
    end
    solutions
end
function scanStark(H::MoleculeHamiltonian, B_field::T, E_field::Vector{<:Real}; Intensity = 1, nev = 36) where T<:Number
    J2Hz = 1.509190311676e33 
    solutions = Vector{sol}(undef, length(E_field))
    for (i, E_i) in enumerate(E_field)
        solutions[i] = diagonalize(H, B_field , E_i, Intensity= 1, nev = nev)
    end
    solutions
end
end

######## Plotting Utils ###############
module plotting
export plotZeemanMap, plotStarkMap, plotIntensityScan
using ..QuantumWrapper.State_Toolbox
using ..QuantumWrapper.AM_Toolbox
import ..Hamiltonian: MoleculeHamiltonian
import ..solve: sol
include("diatomic/plotting.jl")
end

module calculate
export findState, findTransition, transition_dipole_moment, electric_moment, magnetic_moment, diabaticRamp, findMaxOverlap, findAvoidedCrossing
using ..QuantumWrapper.State_Toolbox
import ..QuantumWrapper.AM_Toolbox: Basis, Node, getBasisUC, endNode
import ..Hamiltonian: MoleculeHamiltonian, DipoleMatrix, zeeman_ham, dc

import ..solve: sol
import ..QuantumWrapper.State_Toolbox.State
import ..QuantumWrapper.State_Toolbox.Ket
import ..QuantumWrapper.State_Toolbox.KetName
include("diatomic/calculate.jl")

State(comp::Array{ComplexF64, 2}, QN::Array{Float64, 2}, mol::MoleculeHamiltonian) = State(normalize!(comp), QN, getBasisUC(mol.MolOp.basisTree))
State(QN::Vector{<:Float64}, mol::MoleculeHamiltonian) = State([1.0], [QN], getBasisUC(mol.MolOp.basisTree) )
State(comp::Vector{<:ComplexF64}, mol::MoleculeHamiltonian) = State(comp, getBasisUC(mol.MolOp.basisTree) , getBasisUC(mol.MolOp.basisTree) )
KetName(state::Vector{<:ComplexF64}, mol::MoleculeHamiltonian; QMorder = [5, 3, 2, 1]) = KetName(state, getBasisUC(mol.MolOp.basisTree), QMorder = QMorder )
KetName(state::State, mol::MoleculeHamiltonian; QMorder = [5, 3, 2, 1]) = KetName(Ket(state), getBasisUC(mol.MolOp.basisTree), QMorder = QMorder )
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
end




end