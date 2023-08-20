

using LinearAlgebra

using SparseArrays





function tensor_sum_sparse(A::AbstractSparseMatrix...)
    return blockdiag(A...)
end

⊕(matrices::AbstractSparseMatrix...) = tensor_sum_sparse(matrices...)


raising_operator(j::Node) = raising_operator(j.spin)
x_operator(j::Node)  =  x_operator(j.spin) 
y_operator(j::Node) = y_operator(j.spin)
z_operator(j::Node) = z_operator(j.spin)
function raising_operator(j::Number)
    dim = Int64(2*j + 1)
    J = zeros(Float64, dim, dim)
    
    for (i, mJ) in enumerate(-j:(j-1))
        
        J[i, i + 1] = sqrt(j*(j + 1) - (mJ)*( mJ+ 1))

    end
    return transpose(J)
end

function raising_operator(j::Basis)
    J_o = SparseMatrixCSC[]

    for (i, J) in enumerate(j.shape)
        push!(J_o, sparse(raising_operator(J)))
    end
    ⊕(J_o...)
end

function x_operator(j::Basis)
    J_plus = raising_operator(j)
    J_minus = adjoint(J_plus)
    0.5*(J_plus + J_minus)
    
end

function y_operator(j::Basis)
    J_plus = raising_operator(j)
    J_minus = adjoint(J_plus)
    -0.5*1im*(J_plus - J_minus)
end

function z_operator(j::Basis)

    #⊕(J_o...)
    J_plus = raising_operator(j)
    J_minus = adjoint(J_plus)
    0.5*(J_plus*J_minus - J_minus*J_plus)
end
