module QuantumWrapper

export AM_Toolbox, SpinOperator, State_Toolbox

module AM_Toolbox
include("AM_Toolbox.jl")
export Node, Basis, endNode,  length, SpinBasis, AngularNode, getBasisFC, generateLayers, layerwise_transformation, ⊗, angularMomentum, getUnitary, couple, is_unitary, getBasisUC, getBasis, step_up

end

module SpinOperator
import ..AM_Toolbox: Basis, Node
export ⊕, x_operator, y_operator, z_operator, raising_operator
include("spinOperator.jl")
end

module State_Toolbox
export KetName, State, Ket, +, *, normpath
include("StateUtil.jl")
end
end






