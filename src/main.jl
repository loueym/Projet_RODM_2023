include("building_tree.jl")
include("utilities.jl")

function main(add_cuts::Bool=false)
    if add_cuts
        println("Adding cuts via callback")
    end

    # Pour chaque jeu de données
    # for dataSetName in ["iris", "seeds", "wine", "divorce", "higher_education", "accent"]
    # for dataSetName in ["divorce", "higher_education", "accent"]
    for dataSetName in ["iris"]
        
        print("=== Dataset ", dataSetName)

        # Préparation des données
        include("../data/" * dataSetName * ".txt")

        # Ramener chaque caractéristique sur [0, 1]
        reducedX = Matrix{Float64}(X)
        for j in 1:size(X, 2)
            reducedX[:, j] .-= minimum(X[:, j])
            reducedX[:, j] ./= maximum(X[:, j])
        end

        train, test = train_test_indexes(length(Y))
        X_train = reducedX[train, :]
        Y_train = Y[train]
        X_test = reducedX[test, :]
        Y_test = Y[test]
        classes = unique(Y)

        println(" (train size ", size(X_train, 1), ", test size ", size(X_test, 1), ", ", size(X_train, 2), ", features count: ", size(X_train, 2), ")")
        
        # Temps limite de la méthode de résolution en secondes
        time_limit = 120
        println("Attention : le temps est fixé à $time_limit s.")

        # Pour chaque profondeur considérée
        for D in 2:4

            println("  D = ", D)

            ## 1 - Univarié (séparation sur une seule variable à la fois)
            # Création de l'arbre
            print("    Univarié...  \t")
            T, obj, resolution_time, gap = build_tree(X_train, Y_train, D,  classes, multivariate = false, add_cuts = add_cuts, time_limit = time_limit)

            # Test de la performance de l'arbre
            print(round(resolution_time, digits = 1), "s\t")
            print("gap ", round(gap, digits = 1), "%\t")
            if T != nothing
                print("Erreurs train/test ", prediction_errors(T,X_train,Y_train, classes))
                print("/", prediction_errors(T,X_test,Y_test, classes), "\t")
            end
            println()

            ## 2 - Multivarié
            print("    Multivarié...\t")
            T, obj, resolution_time, gap = build_tree(X_train, Y_train, D, classes, multivariate = true, add_cuts = add_cuts, time_limit = time_limit)
            print(round(resolution_time, digits = 1), "s\t")
            print("gap ", round(gap, digits = 1), "%\t")
            if T != nothing
                print("Erreurs train/test ", prediction_errors(T,X_train,Y_train, classes))
                print("/", prediction_errors(T,X_test,Y_test, classes), "\t")
            end
            println("\n")
        end
    end 
end

# main()
main(true)