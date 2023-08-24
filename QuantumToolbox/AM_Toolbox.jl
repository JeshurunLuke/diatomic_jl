using WignerSymbols
import Base.length
using LinearAlgebra
using SparseArrays
using Base.Threads
using DataStructures


abstract type Basis end
abstract type Node end






struct SpinBasis{T} <: Basis
    QM_Name::String
    shape::Vector{T}
    Coupled::Vector{SpinBasis}
end

struct AngularNode <: Node
    name::String
    spin::SpinBasis
    left::Union{Nothing, AngularNode}
    right::Union{Nothing, AngularNode}
    unitary::Union{Nothing, Matrix{<:Any}, SparseMatrixCSC{<:Any,Int}}
end

function length(a::SpinBasis)
    Int(sum([2*qm + 1 for qm in a.shape]))
end
⊗(a::Vector{Any}, b::Vector{Any})  = vec([vcat(el...) for el in [collect(Iterators.product(a, b))...]])
function tensor_product(vectors::Vector{<:Any}...)
    result = vectors[1]
    for vec in vectors[2:end]
        result = [[r..., v...] for r in result for v in vec]
    end
    return result
end

# Update the tensor function
⊗(vectors::Vector{<:Any}...) = tensor_product(vectors...)






function endNode(root::AngularNode)
    layers = generateLayers(root)
    endNodes = AngularNode[]
    for layer in layers[2:end]
        for node in (layer)
            if node.left == nothing && node.right == nothing
                push!(endNodes, node)
            end
        end
    end
    return endNodes
end
function generateLayers(root::AngularNode)
    layers = [[root]]
    while true
        next_layer = []
        for node in layers[end]
            node.left !== nothing && push!(next_layer, node.left)
            node.right !== nothing && push!(next_layer, node.right)
        end
        if isempty(next_layer)
            break
        end
        push!(layers, next_layer)
    end
    layers
end

#Step up generates the unitary transformation from the base most node (UC basis) to the layer specified
step_up(root::AngularNode) = step_up(root::AngularNode, root.name)
function step_up(root::AngularNode, name::String)
    
    target_size = length(getBasis(root.spin))
    layers = generateLayers(root)
    endLayer = (length(layers) + 1) - (getLayerDepth(root, name))
    transformations = []
    for (depth, layer) in enumerate((reverse(layers)[1:endLayer]))
        for node in reverse(layer)
            # Calculate the size of the identity matrix to tensor with
            unitary_rows = size(node.unitary, 1)
            
            # Ensure that the tensoring is possible
            if target_size % unitary_rows != 0
                throw(DimensionMismatch("Cannot resize unitary of size $unitary_rows to target size $target_size"))
            end
            
            id_size = target_size ÷ unitary_rows
            
            # Construct the transformation for this node using sparse identities
            try
                transformation = kron(
                    spdiagm(0 => ones(Complex{Float64}, id_size)), 
                    node.unitary
                )
                push!(transformations, transformation)
            catch e
                println("Error in constructing transformation at node: $(node.name), depth: $depth")
                rethrow(e)
            end
        end
    end
    return prod(transformations)
end


#Layerwise Transformation generates the unitary from a specified layer to the FC basis
function layerwise_transformation(root::AngularNode, name::String)
    endLayer = getLayerDepth(root, name)
    target_size = length(getBasis(root.spin))
    layers = generateLayers(root)
    transformations = []
    
    for (depth, layer) in enumerate(layers[1:endLayer])
        for node in layer
            # Calculate the size of the identity matrix to tensor with
            unitary_rows = size(node.unitary, 1)
            
            # Ensure that the tensoring is possible
            if target_size % unitary_rows != 0
                throw(DimensionMismatch("Cannot resize unitary of size $unitary_rows to target size $target_size"))
            end
            
            id_size = target_size ÷ unitary_rows
            
            # Construct the transformation for this node using sparse identities
            try
                transformation = kron(
                    spdiagm(0 => ones(Complex{Float64}, id_size)), 
                    node.unitary
                )
                push!(transformations, transformation)
            catch e
                println("Error in constructing transformation at node: $(node.name), depth: $depth")
                rethrow(e)
            end
        end
    end
        
    # Multiply the transformations to get the final one
    return prod(reverse(transformations))
end




function layerwise_transformation(root::AngularNode)
    target_size = length(getBasis(root.spin))
    layers = generateLayers(root)
    transformations = []
    
    for (depth, layer) in enumerate(layers)
        for node in layer
            # Calculate the size of the identity matrix to tensor with
            unitary_rows = size(node.unitary, 1)
            
            # Ensure that the tensoring is possible
            if target_size % unitary_rows != 0
                throw(DimensionMismatch("Cannot resize unitary of size $unitary_rows to target size $target_size"))
            end
            
            id_size = target_size ÷ unitary_rows
            
            # Construct the transformation for this node using sparse identities
            try
                transformation = kron(
                    spdiagm(0 => ones(Complex{Float64}, id_size)), 
                    node.unitary
                )
                push!(transformations, transformation)
            catch e
                println("Error in constructing transformation at node: $(node.name), depth: $depth")
                rethrow(e)
            end
        end
    end
        
    # Multiply the transformations to get the final one
    return prod(reverse(transformations))
end

SpinBasis(QM_Name::String,shape::Vector{T}) where T = SpinBasis{T}(QM_Name, shape, SpinBasis[])
function angularMomentum(QM_Name::String, shape::Vector{T}) where T
    dims=  Int(sum([2*qm + 1 for qm in shape]))
    AngularNode(QM_Name, SpinBasis{T}(QM_Name, shape, SpinBasis[]), nothing, nothing,sparse(1.0I, dims, dims))
end



getBasis(A1::AngularNode, type::String) = getBasis(A1.spin, type)
function getBasisUC(A1::AngularNode)
    vectors_to_tensor = ([getBasis(basis.spin) for basis in endNode(A1)])
    return ⊗(vectors_to_tensor...)
end
getBasisFC(root::AngularNode) = getBasisFC(root, root.name)

function getLayerDepth(root::AngularNode, name::String)
    layers = generateLayers(root)
    endLayer = 1
    for (depth, layer) in enumerate(layers)
        for node in layer
            if node.name== name
                endLayer = depth
            end
        end
    end
    endLayer
end


function getBasisFC(root::AngularNode, name::String)
    endLayer = getLayerDepth(root, name)
    endNodes  = []
    traverse_to_layer!(endNodes, root, endLayer, 1)
    ⊗([getBasis(node, "FC") for node in endNodes]...)
end
function traverse_to_layer!(save_endNodes, start_node::AngularNode, target_layer::Int, current_layer::Int=0)
    if current_layer == target_layer || (start_node.left == nothing && start_node.right == nothing)
        push!(save_endNodes, start_node)  # or any other operation you'd like to perform on the node
        return
    end

    for child in [start_node.left, start_node.right]
        traverse_to_layer!(save_endNodes, child, target_layer, current_layer + 1)
    end
end
function getBasis(a::SpinBasis)
    basis = []
    for basisSpace in a.shape, m in -basisSpace[1]:basisSpace[1]
            push!(basis, vcat([m, vcat(basisSpace...)]...))
    end
    basis
end
function getBasis(root::SpinBasis, type::String)
    if length(root.Coupled) == 2 && (type in ["PC", "FC"])
        getBasis(root.Coupled[1], root.Coupled[2], type)
    elseif length(root.Coupled) == 0
        getBasis(root)
    end

end
function getBasis(B1::SpinBasis, B2::SpinBasis, type::String)
    basis = []

    # Explicit handling for Fully Coupled and Uncoupled basis
    if type == "FC"
        for qm_1 in B1.shape, qm_2 in B2.shape, qm_c in abs(qm_1[1] - qm_2[1]):(qm_1[1] + qm_2[1]),  mq_c in -qm_c:qm_c
            if length(qm_1) > 1&& length(qm_2)  > 1
                push!(basis, [qm_c, mq_c, qm_1[1], qm_2[1], qm_1[2:end]..., qm_2[2:end]...])
            elseif length(qm_1) > 1
                push!(basis, [qm_c, mq_c, qm_1[1], qm_2[1], qm_1[2:end]...])
            elseif  length(qm_2) > 1
                push!(basis, [qm_c, mq_c, qm_1[1], qm_2[1], qm_2[2:end]...])
            else 
                push!(basis, [qm_c, mq_c, qm_1[1], qm_2[1]])

            end
        end
    elseif type == "PC"
        for qm_1 in B1.shape,mq1 in -qm_1[1]:qm_1[1], qm_2 in B2.shape,  mq2 in -qm_2[1]:qm_2[1]
            if length(qm_1) > 1 && length(qm_2) > 1
                push!(basis, [mq1, mq2, qm_1[1], qm_2[1], qm_1[2:end]..., qm_2[2:end]...])
            elseif length(qm_1) > 1
                push!(basis, [mq1, mq2, qm_1[1], qm_2[1], qm_1[2:end]...])
            elseif length(qm_2) > 1
                push!(basis, [mq1, mq2, qm_1[1], qm_2[1], qm_2[2:end]...])
            else 
                push!(basis, [mq1, mq2, qm_1[1], qm_2[1]])
            end

        end
    end
    
    return basis
end


getUnitary(A1::AngularNode, A2::AngularNode) = getUnitary(A1.spin, A2.spin)
function getUnitary(B1::SpinBasis, B2::SpinBasis)
    dim1 = Int(sum([2*qm[1] + 1 for qm in B1.shape]))
    dim2 = Int(sum([2*qm[1] + 1 for qm in B2.shape]))

    #UnitaryMatrix = zeros(ComplexF64, dim1*dim2, dim1*dim2)  # Note the change here

    UCbasis = getBasis(B1, B2, "PC")
    FCbasis = getBasis(B1, B2, "FC")
    
   UnitaryMatrix = basisTransformMatrix( UCbasis, FCbasis)
   #basisTransformMatrixB!(UnitaryMatrix, UCbasis, FCbasis) 
   return UnitaryMatrix
end

struct IndexedBasis
    basis::Array
    index::Int
end
@inline function identify_unique_combinations(basisT_subspace, basisF_subspace)
    return [(bt, bf) for bt in basisT_subspace for bf in basisF_subspace if bt.basis[1] + bt.basis[2] === bf.basis[2]]
end

@inline function compute_clebschgordan_for_combinations(combinations)
    n = length(combinations)
    I_local = Vector{Int}(undef, n)
    J_local = Vector{Int}(undef, n)
    V_local = Vector{ComplexF64}(undef, n)

   @inbounds @simd for i = 1:n
        bt, bf = combinations[i]
        I_local[i] = bt.index
        J_local[i] = bf.index
        V_local[i] = float(WignerSymbols.clebschgordan(bt.basis[3], bt.basis[1], bt.basis[4], bt.basis[2], bf.basis[1], bf.basis[2]))
    end
    return I_local, J_local, V_local
end

function build_subspace_dict(basis)
    dict = Dict()
    for (idx, b) in enumerate(basis)
        subspace = b[3:end]
        indexed_basis = IndexedBasis(b, idx)
        push!(get!(dict, subspace, []), indexed_basis)
    end
    return dict
end

function basisTransformMatrix(basisT, basisF)
    basisT_dict = build_subspace_dict(basisT)
    basisF_dict = build_subspace_dict(basisF)
    subspaces_list = collect(keys(basisT_dict))

    results = Tuple{Vector{Int}, Vector{Int}, Vector{ComplexF64}}[]

    @threads for subspace in subspaces_list
        basisT_subspace = get(basisT_dict, subspace, [])
        basisF_subspace = get(basisF_dict, subspace, [])

        combinations = identify_unique_combinations(basisT_subspace, basisF_subspace)
        push!(results, compute_clebschgordan_for_combinations(combinations))
    end

    I = vcat((result[1] for result in results)...)
    J = vcat((result[2] for result in results)...)
    V = vcat((result[3] for result in results)...)

    return sparse(I, J, V)
end


function couple(basisName::String, A1::AngularNode, A2::AngularNode)
    B1, B2 = A1.spin, A2.spin
    qm_new = []
    for qm_1 in B1.shape,  qm_2 in B2.shape,  qm_c in abs(qm_1[1]-qm_2[1]):(qm_1[1] + qm_2[1])
        push!(qm_new, [qm_c, qm_1, qm_2])
    end
    AngularNode(basisName, SpinBasis(basisName, qm_new, SpinBasis[B1, B2]), A1,A2, getUnitary(B1, B2) )
end


function couple(basisName::String, B1::SpinBasis, B2::SpinBasis)
    qm_new = []
    for qm_1 in B1.shape,  qm_2 in B2.shape,  qm_c in abs(qm_1[1]-qm_2[1]):(qm_1[1] + qm_2[1])
        push!(qm_new, [qm_c, qm_1, qm_2])
    end
    SpinBasis(basisName, qm_new, SpinBasis[B1, B2])
end

function is_unitary(m)
    return isapprox(I, m*adjoint(m))
end 

