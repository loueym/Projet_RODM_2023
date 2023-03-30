# Fonctions utiles permettant la construction de l'arbre avec un callback pour renforcer la formulation 

using Random
using IterTools

# Fonction récupérant les variables en fonction de la formulation choisie pour la séparation
function get_variables_value(model::JuMP.Model, cb_data::CPLEX.CallbackContext)
    a_val = callback_value.(cb_data, model[:a])
    b_val = callback_value.(cb_data, model[:b])
    c_val = callback_value.(cb_data, model[:c])
    u_at_val = callback_value.(cb_data, model[:u_at])
    u_tw_val = callback_value.(cb_data, model[:u_tw])
    return a_val, b_val, c_val, u_at_val, u_tw_val
end

# Fonction permettant de déterminer si c'est l'obtention d'une solution fractionnaire qui a entraîné l'appel d'un callback
function is_fractional_point(cb_data::CPLEX.CallbackContext, context_id::Clong)
    # # Solver-independant way of doing things
    # status = callback_node_status(cb_data, model)
    # if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
    #     # `callback_value(cb_data, x)` is not integer (to some tolerance).
    #     # If, for example, your lazy constraint generator requires an
    #     # integer-feasible primal solution, you can add a `return` here.
    # elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
    #     # `callback_value(cb_data, x)` is integer (to some tolerance).
    #     return
    # else
    #     @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
    #     # `callback_value(cb_data, x)` might be fractional or integer.
    # end

    # context_id  == CPX_CALLBACKCONTEXT_RELAXATION si le  callback est
    # appelé dans le cas ou une nouvelle solution relaxée est obtenue
    if context_id == CPX_CALLBACKCONTEXT_RELAXATION
        return true
    else
        return false
    end
end

# Fonction permettant l'ajout d'inégalité valides au modèle
function add_user_cut(model::JuMP.Model, cb_data::CPLEX.CallbackContext, x::Matrix{Float64}, y::Vector, multivariate::Bool=false, verbose::Bool=false)
    # Récupération des valeurs et des variables
    a, b, c, u_at, u_tw = model[:a], model[:b], model[:c], model[:u_at], model[:u_tw]
    a_val, b_val, c_val, u_at_val, u_tw_val = get_variables_value(model, cb_data)
    # Mettre le méchanisme de sélection de l'inégalité choisie
    J, N = size(a_val)
    I, K = size(u_at_val)[1], size(c_val)[1]
    index_J, index_N, index_K, index_I = shuffle(1:J), shuffle(1:N), shuffle(1:K), shuffle(1:I)
    # Ajout de la coupe 2 si elle a été trouvée
    max_cut, k, i, t = secondCut(index_N, index_K, index_I, y, c_val, u_at_val, u_tw_val)
    if max_cut > 1
        for i in 1:I
            if y[i] != k
                cut = @build_constraint(c[k, t] + u_tw[i, t] + u_at[i, t*2] + u_at[i, t*2+1] <= 1)
                MOI.submit(model, MOI.UserCut(cb_data), cut)
            end
        end
    end
    if verbose
        println("Adding valid inequality : $(cut)")
    end
    # if !multivariate
    #     index_I1_I2 = collect(subsets(index_I, 2))
    #     # # Ajout de la coupe (1) si elle a été trouvée
    #     # max_cut, i1, i2, t = firstCut(index_N, index_J, index_I1_I2, x, u_at_val)
    #     # if max_cut > 1
    #     #     cut = @build_constraint(u_at[i1, t*2] + u_at[i2, t*2+1] <= 1)
    #     #     MOI.submit(model, MOI.UserCut(cb_data), cut)
    #     # end
    #     # if verbose
    #     #     println("Adding valid inequality : $(cut)")
    #     # end
    #     # Ajout de la coupe (3) si elle a été trouvée
    #     max_cut, i1, i2, t = thirdCut(index_N, K, J, index_I1_I2, x, c_val, u_at_val, a_val)
    #     if max_cut > 2
    #         cut = @build_constraint(2 * sum(c[k, t] for k in 1:K) + u_at[i1, t*2+1] + u_at[i2, t*2] + sum(a[j, t] for j in 1:J if x[i1, j] <= x[i2, j]) <= 2)
    #         MOI.submit(model, MOI.UserCut(cb_data), cut)
    #         if verbose
    #             println("Adding valid inequality : $(cut)")
    #         end
    #     end
    # end
end

# Fonction qui retourne la première inégalité (1) violée ou 0 sinon
function firstCut(index_N::Vector{Int}, index_J::Vector{Int}, index_I1_I2::Vector{Vector{Int}}, x::Matrix{Float64}, u_at_star)
    for t in index_N
        for j in index_J
            for (i1, i2) in index_I1_I2
                if x[i1, j] >= x[i2, j]
                    current_cut = u_at_star[i1, t*2] + u_at_star[i2, t*2+1]
                    if current_cut > 1 + 1e-5
                        return current_cut, i1, i2, t
                    end
                else
                    current_cut = u_at_star[i1, t*2+1] + u_at_star[i2, t*2]
                    if current_cut > 1 + 1e-5
                        return current_cut, i2, i1, t
                    end
                end
            end
        end
    end
    return 0, 0, 0, 0
end

# Fonction qui retroune la première inégalité (2) violée ou 0 sinon
function secondCut(index_N::Vector{Int}, index_K::Vector{Int}, index_I::Vector{Int}, y::Vector, c_star, u_at_star, u_tw_star)
    for i in index_I
        for k in index_K
            for t in index_N
                if y[i] != k
                    current_cut = c_star[k, t] + u_tw_star[i, t] + u_at_star[i, t*2] + u_at_star[i, t*2+1]
                    if current_cut > 1 + 1e-5
                        return current_cut, k, i, t
                    end
                end
            end
        end
    end
    # Testing mean deviation
    return 0, 0, 0, 0
end

# Fonction qui retroune la première inégalité (3) violée ou 0 sino
function thirdCut(index_N::Vector{Int}, K::Int, J::Int, index_I1_I2::Vector{Vector{Int}}, x::Matrix{Float64}, c_star, u_at_star, a_star)
    for t in index_N
        for (i1, i2) in index_I1_I2
            a_sep = [a_star[j, t] for j in 1:J if x[i1, j] <= x[i2, j]]
            if length(a_sep) > 0
                current_cut = 2 * sum(c_star[k, t] for k in 1:K) + u_at_star[i1, t*2+1] + u_at_star[i2, t*2] + sum(a_sep)
                if current_cut > 2 + 1e-5
                    return current_cut, i1, i2, t
                end
            end
        end
    end  
    return 0, 0, 0, 0  
end

# Fonction pour avoir t et ses ancêtres
function get_t_and_ancestors(t::Int)
    t_and_ancestors = [t]
    while t > 1
        t = floor(Int, t/2)
        push!(t_and_ancestors, t)
    end
    return t_and_ancestors
end
