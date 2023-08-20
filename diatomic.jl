module diatomic_jl


include("QuantumToolbox/QuantumWrapper.jl")
using .QuantumWrapper
import .QuantumWrapper: AM_Toolbox
export Molecule, Hamiltonian, solve, plotting, calculate, AM_Toolbox, getBasisUC

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
export generateHamiltonian, OpticalBeam, MoleculeHamiltonian
import ..Molecule: moleculeProperties
include("diatomic/Hamiltonian.jl")
using ..diatomic_jl
import ..diatomic_jl: getBasisUC
getBasisUC(mol::MoleculeHamiltonian) = AM_Toolbox.getBasisUC(mol.MolOp.basisTree)
end

##### Used to solve the Hamiltonian ################################
module solve
export scanZeeman, scanStark, diagonalize, sol
import ..Hamiltonian: MoleculeHamiltonian
using LinearAlgebra
using Arpack

struct sol
    B_field
    E_field
    vec
    val
end

function diagonalize(H::MoleculeHamiltonian, B_field, E_field; Intensity = 1,  nev = 36)
    J2Hz = 1.509190311676e33 
    
    vals, vecs = eigen(Hermitian(Matrix(H.Hhfs .+ H.Hzeem*B_field .+ E_field*H.Hdc + sum(Intensity .* H.Hac))))
    sorted_indices = sortperm(vals)
    vals = vals[sorted_indices]
    vecs = vecs[:, sorted_indices]
    sol(B_field, E_field, vecs, vals*J2Hz)
end
function scanZeeman(H::MoleculeHamiltonian, B_field::Vector{<:Real}, E_field::T) where T<:Number
    J2Hz = 1.509190311676e33 
    solutions = Vector{sol}(undef, length(B_field))
    for (i, B_i) in enumerate(B_field)
        solutions[i] = diagonalize(H, B_i , E_field)
    end
    solutions
end
function scanStark(H::MoleculeHamiltonian, B_field::T, E_field::Vector{<:Real}) where T<:Number
    J2Hz = 1.509190311676e33 
    solutions = Vector{sol}(undef, length(E_field))
    for (i, E_i) in enumerate(E_field)
        solutions[i] = diagonalize(H, B_field , E_i)
    end
    solutions
end
end

######## Plotting Utils ###############
module plotting
export plotZeemanMap, plotStarkMap
using ..QuantumWrapper.State_Toolbox
using ..QuantumWrapper.AM_Toolbox
import ..Hamiltonian: MoleculeHamiltonian
import ..solve: sol
include("diatomic/plotting.jl")
end

module calculate
export findState, findTransition, transition_dipole_moment
using ..QuantumWrapper.State_Toolbox
import ..QuantumWrapper.AM_Toolbox: Basis, Node, getBasisUC
import ..Hamiltonian: MoleculeHamiltonian
import ..solve: sol
import ..QuantumWrapper.State_Toolbox.State
include("diatomic/calculate.jl")

State(comp::Array{ComplexF64, 2}, QN::Array{Float64, 2}, mol::MoleculeHamiltonian) = State(normalize!(comp), QN, getBasisUC(mol.MolOp.basisTree))
State(QN::Vector{<:Any}, mol::MoleculeHamiltonian) = State([1.0], QN, getBasisUC(mol.MolOp.basisTree) )
end




end