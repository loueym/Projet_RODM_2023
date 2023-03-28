include("struct/tree.jl")

"""
Construit un arbre de décision par résolution d'un programme linéaire en nombres entiers

Entrées :
- x : caractéristiques des données d'entraînement
- y : classe des données d'entraînement
- D : Nombre maximal de séparations d'une branche (profondeur de l'arbre - 1)
- multivariate (optionnel): vrai si les séparations sont multivariées; faux si elles sont univariées (faux par défaut)
- mu (optionnel, utilisé en multivarié): distance minimale à gauche d'une séparation où aucune donnée ne peut se trouver (i.e., pour la séparation ax <= b, il n'y aura aucune donnée dans ]b - ax - mu, b - ax[) (10^-4 par défaut)
- time_limits (optionnel) : temps maximal de résolution (-1 si le temps n'est pas limité) (-1 par défaut)
- classes : labels des classes figurant dans le dataset
"""
function build_tree(x::Matrix{Float64}, y::Vector, D::Int64, classes; multivariate::Bool=false, add_cuts::Bool=false, time_limit::Int64 = -1, mu::Float64=10^(-4))
    
    dataCount = length(y) # Nombre de données d'entraînement
    featuresCount = length(x[1, :]) # Nombre de caractéristiques
    classCount = length(classes) # Nombre de classes différentes
    sepCount = 2^D - 1 # Nombre de séparations de l'arbre
    leavesCount = 2^D # Nombre de feuilles de l'arbre

    m = Model(CPLEX.Optimizer) 
    MOI.set(m, MOI.NumberOfThreads(), 1)
    set_silent(m)

    if time_limit!=-1
        set_time_limit_sec(m, time_limit)
    end

    # Plus petite différence entre deux données pour une caractéristique
    mu_min = 1.0 
    # Plus grande différence entre deux données pour une caractéristique
    mu_max = 0.0
    
    if !multivariate # calcul des constantes mu_min, mu_max et du vecteur mu

        # mu_vect[j] est la plus petite différence (>0) entre deux données, pour la caractéristiques j
        mu_vect = ones(Float64, featuresCount)
        for j in 1:featuresCount
            for i1 in 1:dataCount
                for i2 in (i1+1):dataCount
                    if x[i1, j] - x[i2, j] != 0
                        mu_vect[j] = min(mu_vect[j], abs(x[i1, j] - x[i2, j]))
                    end
                end
            end
            mu_min = min(mu_min, mu_vect[j])
            mu_max = max(mu_max, mu_vect[j])
        end
    end

    ## Déclaration des variables
    if multivariate
        @variable(m, a[1:featuresCount, 1:sepCount], base_name="a_{j, t}")
        @variable(m, a_h[1:featuresCount, 1:sepCount], base_name="â_{j, t}")
        @variable(m, s[1:featuresCount, 1:sepCount], Bin, base_name="s_{j, t}")
        @variable(m, d[1:sepCount], Bin, base_name="d_t")
    else
        @variable(m, a[1:featuresCount, 1:sepCount], Bin, base_name="a")
    end 
    @variable(m, b[1:sepCount], base_name="b_t")
    @variable(m, c[1:classCount, 1:(sepCount+leavesCount)], Bin, base_name = "c_{k, t}")
    @variable(m, u_at[1:dataCount, 1:(sepCount+leavesCount)], Bin, base_name = "u^i_{a(t), t}")
    @variable(m, u_tw[1:dataCount, 1:(sepCount+leavesCount)], Bin, base_name = "u^i_{t, w}")

    ## Déclaration des contraintes

    # contraintes définissant la structure de l'arbre
    if multivariate
        @constraint(m, [t in 1:sepCount], d[t] + sum(c[k, t] for k in 1:classCount) == 1) # on s'assure que le noeud applique une règle de branchement OU attribue une classe
        @constraint(m, [t in 1:sepCount], b[t] <= d[t]) # b doit être nul si il n'y a pas de branchement 
        @constraint(m, [t in 1:sepCount], b[t] >= -d[t]) # b doit être nul si il n'y a pas de branchement 
        @constraint(m, [t in 1:sepCount], sum(a_h[j, t] for j in 1:featuresCount) <= d[t]) # on borne la norme du vecteur a
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], a[j, t] <= a_h[j, t]) # définition de â borne sup de la valeur absolu de a
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], a[j, t] >= -a_h[j, t]) # définition de â borne sup de la valeur absolu de a
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], a[j, t] <= s[j, t]) # définition de s, non nul ssi a non nul
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], a[j, t] >= -s[j, t]) # définition de s, non nul ssi a non nul
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], s[j, t] <= d[t]) # définition de d, non nul si il existe un coef non nul
        @constraint(m, [t in 1:sepCount], sum(s[j, t] for j in 1:featuresCount) >= d[t]) # définition de d, non nul si il existe un coef non nul
        @constraint(m, [t in 2:sepCount], d[t] <= d[t ÷ 2]) # on s'assure que si un noeud de branchement n'applique pas de règle de branchement, ses fils non plus
    else
        @constraint(m, [t in 1:sepCount], sum(a[j, t] for j in 1:featuresCount) + sum(c[k, t] for k in 1:classCount) == 1) # on s'assure que le noeud applique une règle de branchement OU attribue une classe
        @constraint(m, [t in 1:sepCount], b[t] <= sum(a[j, t] for j in 1:featuresCount)) # b doit être nul si il n'y a pas de branchement 
        @constraint(m, [t in 1:sepCount], b[t] >= 0) # b doit être positif
    end
    @constraint(m, [t in (sepCount+1):(sepCount+leavesCount)], sum(c[k, t] for k in 1:classCount) == 1) # on s'assure qu'on attribue une classe par feuille

    # contraintes de conservation du flot et contraintes de capacité
    @constraint(m, [i in 1:dataCount, t in 1:sepCount], u_at[i, t] == u_at[i, t*2] + u_at[i, t*2+1] + u_tw[i, t]) # conservation du flot dans les noeuds de branchement
    @constraint(m, [i in 1:dataCount, t in (sepCount+1):(sepCount+leavesCount)], u_at[i, t] == u_tw[i, t]) # conservation du flot dans les feuilles
    @constraint(m, [i in 1:dataCount, t in 1:(sepCount+leavesCount)], u_tw[i, t] <= c[findfirst(classes .== y[i]), t]) # contrainte de capacité qui impose le flot a etre nul si la classe de la feuille n'est pas la bonne
    if multivariate
        @constraint(m, [i in 1:dataCount, t in 1:sepCount], sum(a[j, t]*x[i, j] for j in 1:featuresCount) + mu <= b[t] + (2+mu)*(1-u_at[i, t*2])) # contrainte de capacité controlant le passage dans le noeud fils gauche
        @constraint(m, [i in 1:dataCount, t in 1:sepCount], sum(a[j, t]*x[i, j] for j in 1:featuresCount) >= b[t] - 2*(1-u_at[i, t*2 + 1])) # contrainte de capacité controlant le passage dans le noeud fils droit
        @constraint(m, [i in 1:dataCount, t in 1:sepCount], u_at[i, t*2+1] <= d[t]) # contrainte de capacité empechant les données de passer dans le fils droit d'un noeud n'appliquant pas de règle de branchement
    else
        @constraint(m, [i in 1:dataCount, t in 1:sepCount], sum(a[j, t]*(x[i, j]+mu_vect[j]-mu_min) for j in 1:featuresCount) + mu_min <= b[t] + (1+mu_max)*(1-u_at[i, t*2])) # contrainte de capacité controlant le passage dans le noeud fils gauche
        @constraint(m, [i in 1:dataCount, t in 1:sepCount], sum(a[j, t]*x[i, j] for j in 1:featuresCount) >= b[t] - (1-u_at[i, t*2 + 1])) # contrainte de capacité controlant le passage dans le noeud fils droit
        @constraint(m, [i in 1:dataCount, t in 1:sepCount], u_at[i, t*2+1] <= sum(a[j, t] for j in 1:featuresCount)) # contrainte de capacité empechant les données de passer dans le fils droit d'un noeud n'appliquant pas de règle de branchement
        @constraint(m, [i in 1:dataCount, t in 1:sepCount], u_at[i, t*2] <= sum(a[j, t] for j in 1:featuresCount)) # contrainte de capacité empechant les données de passer dans le fils gauche d'un noeud n'appliquant pas de règle de branchement
    end

    ## Déclaration de l'objectif
    if multivariate
        @objective(m, Max,  sum(u_at[i, 1] for i in 1:dataCount))
    else
        @objective(m, Max, sum(u_at[i, 1] for i in 1:dataCount))
    end

    classif = @expression(m, sum(u_at[i, 1] for i in 1:dataCount))

    
    function callback_function(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if isIntegerPoint(cb_data, context_id)
            callback_called = true
            CPLEX.load_callback_variable_primal(cb_data, context_id)

            if multivariate
                a_star = callback_value.(cb_data, s)
            else
                a_star = callback_value.(cb_data, a)
            end
            c_star = callback_value.(cb_data, c)
            u_at_star = callback_value.(cb_data, u_at)
            u_tw_star = callback_value.(cb_data, u_tw)

            max_cut, i1, i2, t = firstCut(sepCount, featuresCount, dataCount, x, u_at_star)
            if max_cut > 1
                cut = @build_constraint(u_at[i1, t*2] + u_at[i2, t*2+1] <= 1)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cut)
            end

            max_cut, k, i, t = secondCut(sepCount, classCount, dataCount, y, c_star, u_at_star, u_tw_star)
            if max_cut > 1
                cut = @build_constraint(c[k, t] + u_tw[i, t] + u_at[i, t*2] + u_at[i, t*2+1] <= 1)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cut)
            end

            max_cut, i1, i2, t = thirdCut(sepCount, classCount, featuresCount, dataCount, x, c_star, u_at_star, a_star)
            if max_cut > 2
                cut = @build_constraint(2 * sum(c[k, t] for k in 1:classCount) + u_at[i1, t*2+1] + u_at[i2, t*2] + sum(a[j, t] for j in 1:featuresCount if x[i1, j] <= x[i2, j]) <= 2)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cut)
            end
        end
    end

    if add_cuts
        MOI.set(m, CPLEX.CallbackFunction(), callback_function)
    end

    starting_time = time()
    optimize!(m)
    resolution_time = time() - starting_time
    
    gap = -1.0

    # Arbre obtenu (vide si le solveur n'a trouvé aucune solution)
    T = nothing
    
    # Si une solution a été trouvée
    if primal_status(m) == MOI.FEASIBLE_POINT

        # class[t] : classe prédite par le sommet t
        class = Vector{Int64}(undef, sepCount+leavesCount)
        for t in 1:(sepCount+leavesCount)
            k = argmax(value.(c[:, t]))
            if value.(c[k, t])  >=  1.0 - 10^-4
                class[t] = k
            else
                class[t] = -1
            end
        end
        
        # Si une solution optimale a été trouvée
        if termination_status(m) == MOI.OPTIMAL
            gap = 0
        else
            # Calcul du gap relatif entre l'objectif de la meilleure solution entière et la borne continue en fin de résolution
            objective = JuMP.objective_value(m)
            bound = JuMP.objective_bound(m)
            gap = 100.0 * abs(objective - bound) / (objective + 10^-4) # +10^-4 permet d'éviter de diviser par 0
        end   
        
        # Construction d'une variable de type Tree dans laquelle chaque séparation est recentrée
        if multivariate
            T = Tree(D, class, round.(Int, value.(u_at)), round.(Int, value.(s)), x)
        else
            T = Tree(D, value.(a), class, round.(Int, value.(u_at)), x)
        end
    end   

    return T, objective_value(m), resolution_time, gap
end


function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)
    # S’il n’y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end

function firstCut(N::Int, J::Int, I::Int, x::Matrix{Float64}, u_at_star)
    max_cut = 0
    best_i1, best_i2, best_t = 0, 0, 0, 0
    for t in 1:N
        for j in 1:J
            for i1 in 1:I
                for i2 in i1+1:I
                    if x[i1, j] >= x[i2, j]
                        current_cut = u_at_star[i1, t*2] + u_at_star[i2, t*2+1]
                        if current_cut > max_cut
                            max_cut = current_cut
                            best_i1, best_i2, best_t = i1, i2, t
                        end
                    else
                        current_cut = u_at_star[i1, t*2+1] + u_at_star[i2, t*2]
                        if current_cut > max_cut
                            max_cut = current_cut
                            best_i1, best_i2, best_t = i2, i1, t
                        end
                    end
                end
            end
        end
    end
    return max_cut, best_i1, best_i2, best_t
end

function secondCut(N::Int, K::Int, I::Int, y::Vector, c_star, u_at_star, u_tw_star)
    max_cut = 0
    best_k, best_i, best_t = 0, 0, 0
    for t in 1:N
        for k in 1:K
            for i in 1:I
                if y[i] != k
                    current_cut = c_star[k, t] + u_tw_star[i, t] + u_at_star[i, t*2] + u_at_star[i, t*2+1]
                    if current_cut > max_cut
                        max_cut = current_cut
                        best_k, best_i, best_t = k, i, t
                    end
                end
            end
        end
    end
    return max_cut, best_k, best_i, best_t
end


function thirdCut(N::Int, K::Int, J::Int, I::Int, x::Matrix{Float64}, c_star, u_at_star, a_star)
    max_cut = 0
    best_i1, best_i2, best_t = 0, 0, 0
    for t in 1:N
        for i1 in 1:I
            for i2 in 1:I
                if i1 != i2
                    a_sep = [a_star[j, t] for j in 1:J if x[i1, j] <= x[i2, j]]
                    if length(a_sep) > 0
                        current_cut = 2 * sum(c_star[k, t] for k in 1:K) + u_at_star[i1, t*2+1] + u_at_star[i2, t*2] + sum(a_sep)
                        if current_cut > max_cut
                            max_cut = current_cut
                            best_i1, best_i2, best_t = i1, i2, t
                        end
                    end
                end
            end
        end
    end  
    return max_cut, best_i1, best_i2, best_t  
end

"""
FONCTION SIMILAIRE A LA PRECEDENTE UTILISEE UNIQUEMENT SI VOUS FAITES DES REGROUPEMENTS 

Construit un arbre de décision par résolution d'un programme linéaire en nombres entiers

Entrées :
- clusters : partition des données d'entraînement (chaque cluster contient des données de même classe)
- D : Nombre maximal de séparations d'une branche (profondeur de l'arbre - 1)
- multivariate (optionnel): vrai si les séparations sont multivariées; faux si elles sont univariées (faux par défaut)
- mu (optionnel, utilisé en multivarié): distance minimale à gauche d'une séparation où aucune donnée ne peut se trouver (i.e., pour la séparation ax <= b, il n'y aura aucune donnée dans ]b - ax - mu, b - ax[) (10^-4 par défaut)
- time_limits (optionnel) : temps maximal de résolution (-1 si le temps n'est pas limité) (-1 par défaut)
"""
function build_tree(clusters::Vector{Cluster}, D::Int64, classes;multivariate::Bool=false, time_limit::Int64 = -1, mu::Float64=10^(-4))
    
    clusterCount = length(clusters) # Nombre de données d'entraînement
    featuresCount = length(clusters[1].lBounds) # Nombre de caractéristiques
    classCount = length(classes) # Nombre de classes différentes
    sepCount = 2^D - 1 # Nombre de séparations de l'arbre
    leavesCount = 2^D # Nombre de feuilles de l'arbre
    
    m = Model(CPLEX.Optimizer) 

    set_silent(m) # Masque les sorties du solveur

    if time_limit!=-1
        set_time_limit_sec(m, time_limit)
    end

    # Plus petite différence entre deux données pour une caractéristique
    mu_min = 1.0 
    # Plus grande différence entre deux données pour une caractéristique
    mu_max = 0.0
    
    if !multivariate # calcul des constantes mu_min, mu_max et du vecteur mu
        mu_vect = ones(Float64, featuresCount)
        for j in 1:featuresCount
            for i1 in 1:clusterCount
                for i2 in (i1+1):clusterCount

                    v1 = clusters[i1].lBounds[j] - clusters[i2].uBounds[j]
                    v2 = clusters[i2].lBounds[j] - clusters[i1].uBounds[j]

                    # Si les clusters n'ont pas des intervalles pour la caractéristique j qui s'intersectent
                    if v1 > 0 || v2 > 0
                        vMin = min(abs(v1), abs(v2))
                        mu_vect[j] = min(mu_vect[j], vMin)
                    end
                end
            end
            mu_min = min(mu_min, mu_vect[j])
            mu_max = max(mu_max, mu_vect[j])
        end
    end

    ## Déclaraction des variables
    if multivariate
        @variable(m, a[1:featuresCount, 1:sepCount], base_name="a_{j, t}")
        @variable(m, a_h[1:featuresCount, 1:sepCount], base_name="â_{j, t}")
        @variable(m, s[1:featuresCount, 1:sepCount], Bin, base_name="s_{j, t}")
        @variable(m, d[1:sepCount], Bin, base_name="d_t")
    else
        @variable(m, a[1:featuresCount, 1:sepCount], Bin, base_name="a")
    end 
    @variable(m, b[1:sepCount], base_name="b_t")
    @variable(m, c[1:classCount, 1:(sepCount+leavesCount)], Bin, base_name = "c_{k, t}")
    @variable(m, u_at[1:clusterCount, 1:(sepCount+leavesCount)], Bin, base_name = "u^i_{a(t), t}")
    @variable(m, u_tw[1:clusterCount, 1:(sepCount+leavesCount)], Bin, base_name = "u^i_{t, w}")

    ## Déclaration des contraintes
    
    # Contraintes définissant la structure de l'arbre
    if multivariate
        @constraint(m, [t in 1:sepCount], d[t] + sum(c[k, t] for k in 1:classCount) == 1) # on s'assure que le noeud applique une règle de branchement OU attribue une classe
        @constraint(m, [t in 1:sepCount], b[t] <= d[t]) # b doit être nul si il n'y a pas de branchement 
        @constraint(m, [t in 1:sepCount], b[t] >= -d[t]) # b doit être nul si il n'y a pas de branchement 
        @constraint(m, [t in 1:sepCount], sum(a_h[j, t] for j in 1:featuresCount) <= d[t]) # on borne la norme du vecteur a
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], a[j, t] <= a_h[j, t]) # définition de â borne sup de la valeur absolu de a
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], a[j, t] >= -a_h[j, t]) # définition de â borne sup de la valeur absolu de a
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], a[j, t] <= s[j, t]) # définition de s, non nul ssi a non nul
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], a[j, t] >= -s[j, t]) # définition de s, non nul ssi a non nul
        @constraint(m, [t in 1:sepCount, j in 1:featuresCount], s[j, t] <= d[t]) # définition de d, non nul si il existe un coef non nul
        @constraint(m, [t in 1:sepCount], sum(s[j, t] for j in 1:featuresCount) >= d[t]) # définition de d, non nul si il existe un coef non nul
        @constraint(m, [t in 2:sepCount], d[t] <= d[t ÷ 2]) # on s'assure que si un noeud de branchement n'applique pas de règle de branchement, ses fils non plus
    else
        @constraint(m, [t in 1:sepCount], sum(a[j, t] for j in 1:featuresCount) + sum(c[k, t] for k in 1:classCount) == 1) # on s'assure que le noeud applique une règle de branchement OU attribue une classe
        @constraint(m, [t in 1:sepCount], b[t] <= sum(a[j, t] for j in 1:featuresCount)) # b doit être nul si il n'y a pas de branchement 
        @constraint(m, [t in 1:sepCount], b[t] >= 0) # b doit être positif
    end
    @constraint(m, [t in (sepCount+1):(sepCount+leavesCount)], sum(c[k, t] for k in 1:classCount) == 1) # on s'assure qu'on attribue une classe par feuille

    # contraintes de conservation du flot et contraintes de capacité
    @constraint(m, [i in 1:clusterCount, t in 1:sepCount], u_at[i, t] == u_at[i, t*2] + u_at[i, t*2+1] + u_tw[i, t]) # conservation du flot dans les noeuds de branchement
    @constraint(m, [i in 1:clusterCount, t in (sepCount+1):(sepCount+leavesCount)], u_at[i, t] == u_tw[i, t]) # conservation du flot dans les feuilles
    @constraint(m, [i in 1:clusterCount, t in 1:(sepCount+leavesCount)], u_tw[i, t] <= c[findfirst(classes .== clusters[i].class), t]) # contrainte de capacité qui impose le flot a etre nul si la classe de la feuille n'est pas la bonne
    if multivariate
        
        @constraint(m, [i in 1:clusterCount, t in 1:sepCount, dataId in 1:size(clusters[i].x, 1)], sum(a[j, t]*clusters[i].x[dataId, j] for j in 1:featuresCount) + mu <= b[t] + (2+mu)*(1-u_at[i, t*2])) # contrainte de capacité controlant le passage dans le noeud fils gauche
        @constraint(m, [i in 1:clusterCount, t in 1:sepCount, dataId in 1:size(clusters[i].x, 1)], sum(a[j, t]*clusters[i].x[dataId, j] for j in 1:featuresCount) >= b[t] - 2*(1-u_at[i, t*2 + 1])) # contrainte de capacité controlant le passage dans le noeud fils droit
        @constraint(m, [i in 1:clusterCount, t in 1:sepCount], u_at[i, t*2+1] <= d[t]) # contrainte de capacité empechant les données de passer dans le fils droit d'un noeud n'appliquant pas de règle de branchement
    else
        @constraint(m, [i in 1:clusterCount, t in 1:sepCount], sum(a[j, t]*(clusters[i].uBounds[j]+mu_vect[j]-mu_min) for j in 1:featuresCount) + mu_min <= b[t] + (1+mu_max)*(1-u_at[i, t*2])) # contrainte de capacité controlant le passage dans le noeud fils gauche
        @constraint(m, [i in 1:clusterCount, t in 1:sepCount], sum(a[j, t]*clusters[i].lBounds[j] for j in 1:featuresCount) >= b[t] - (1-u_at[i, t*2 + 1])) # contrainte de capacité controlant le passage dans le noeud fils droit
        @constraint(m, [i in 1:clusterCount, t in 1:sepCount], u_at[i, t*2+1] <= sum(a[j, t] for j in 1:featuresCount)) # contrainte de capacité empechant les données de passer dans le fils droit d'un noeud n'appliquant pas de règle de branchement
    end

    ## Déclaration de l'objectif
    if multivariate
        @objective(m, Max, sum(length(clusters[i].dataIds) * u_at[i, 1] for i in 1:clusterCount))
    else
        @objective(m, Max, sum(length(clusters[i].dataIds) * u_at[i, 1] for i in 1:clusterCount))
    end

    starting_time = time()
    optimize!(m)
    resolution_time = time() - starting_time

    # class[t] : classe prédite par le sommet t
    class = Vector{Int64}(undef, sepCount+leavesCount)
    for t in 1:(sepCount+leavesCount)
        k = argmax(value.(c[:, t]))
        if value.(c[k, t]) >= 1.0 - 10^-4
            class[t] = k
        else
            class[t] = -1
        end
    end
        
    gap = -1.0

    # Arbre obtenu (vide si le solveur n'a trouvé aucune solution)
    T = nothing
    
    # Si une solution a été trouvée
    if primal_status(m) == MOI.FEASIBLE_POINT

        # Si une solution optimale a été trouvée
        if termination_status(m) == MOI.OPTIMAL
            gap = 0
        else
            # Calcul du gap relatif entre l'objectif de la meilleure solution entière et la borne continue en fin de résolution
            objective = JuMP.objective_value(m)
            bound = JuMP.objective_bound(m)
            gap = 100.0 * abs(objective - bound) / (objective + 10^-4) # +10^-4 permet d'éviter de diviser par 0
        end   
        
        # Construction d'une variable de type Tree dans laquelle chaque séparation est recentrée
        if multivariate
            T = Tree(D, class, round.(Int, value.(u_at)), round.(Int, value.(s)), clusters)
        else
            T = Tree(D, value.(a), class, round.(Int, value.(u_at)), clusters)
        end
    end   

    return T, objective_value(m), resolution_time, gap
end

