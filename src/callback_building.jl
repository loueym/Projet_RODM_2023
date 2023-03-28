# Fonctions utiles permettant la construction de l'arbre avec un callback pour renforcer la formulation 

# Fonction récupérant les variables en fonction de la formulation choisie pour la séparation
function get_variables_value(model::JuMP.Model, cb_data::CPLEX.CallbackContext, multivariate::Bool=false,)
    a_val = callback_value.(cb_data, model[:a])
    b_val = callback_value.(cb_data, model[:b])
    c_val = callback_value.(cb_data, model[:c])
    u_at_val = callback_value.(cb_data, model[:u_at])
    u_tw_val = callback_value.(cb_data, model[:u_tw])
    if multivariate
        a_h_val = callback_value.(cb_data, model[:a_h])
        s_val = callback_value.(cb_data, model[:s])
        d_val = callback_value.(cb_data, model[:d])
        return a_val, a_h_val, s_val, d_val, b_val, c_val, u_at_val, u_tw_val
    else
        return a_val, nothing, nothing, nothing, b_val, c_val, u_at_val, u_tw_val
    end    
end

# Fonction permettant de déterminer si c'est l'obtention d'une solution entière qui a entraîné l'appel d'un callback
function is_integer_point(cb_data::CPLEX.CallbackContext, context_id::Clong)

    # context_id  == CPX_CALLBACKCONTEXT_CANDIDATE si le  callback est
    # appelé dans un des deux cas suivants :
    # cas 1 - une solution entière a été obtenue; ou
    # cas 2 - une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end

    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)

    # S'il n'y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end

# Fonction permettant l'ajout d'inégalité valides au modèle
function add_lazy_constraint(model::JuMP.Model, cb_data::CPLEX.CallbackContext, multivariate::Bool=false, verbose::Bool=false)
    # Récupération des valeurs des variables
    a_val, a_h_val, s_val, d_val, b_val, c_val, u_at_val, u_tw_val = get_variables_value(model, cb_data, multivariate)
    # Mettre le méchanisme de sélection de la contrainte choisie
    # Compléter par la formulation de la contrainte 
    cstr = @build_constraint()
    MOI.submit(model, MOI.LazyConstraint(cb_data), cstr)
    if verbose
        println("Adding new constraint : $(cstr)")
    end
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
function add_user_cut(model::JuMP.Model, cb_data::CPLEX.CallbackContext, multivariate::Bool=false, verbose::Bool=false)
    # Récupération des valeurs des variables
    a_val, a_h_val, s_val, d_val, b_val, c_val, u_at_val, u_tw_val = get_variables_value(model, cb_data, multivariate)
    # Mettre le méchanisme de sélection de l'inégalité choisie
    # Compléter par la formulation de l'inégalité valide
    cstr = @build_constraint()
    MOI.submit(model, MOI.UserCut(cb_data), cstr)
    if verbose
        println("Adding valid inequality : $(cstr)")
    end
end
