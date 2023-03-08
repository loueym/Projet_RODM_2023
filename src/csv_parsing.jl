# Transforming csv data file to suitable txt format

using CSV 

function transform_csv(file_name::String, y_in_first_column::Bool=false)
    file_path = string(dirname(@__DIR__), "\\data\\", file_name)
    # Reading .csv file
    csv_reader = CSV.File(file_path)
    n = length(csv_reader)
    p = length(csv_reader.names) - 1
    X = zeros(n, p)
    Y = zeros(n)
    println("CSV file $(file_name) parsed (n = $n, p = $p)")
    println("Column names : $(csv_reader.names)")
    for (i, row) in enumerate(csv_reader)
        if y_in_first_column
            Y[i] = row[1]
            for j in 2:p+1
                X[i, j-1] = row[j]
            end
        else
            for j in 1:p
                X[i, j] = row[j]
            end
            Y[i] = row[p + 1]
        end
    end
    # Extracting name 
    file_name = splitext(file_name)[1]
    # Writing .txt file
    open(string(dirname(@__DIR__), "\\data\\", file_name, ".txt"), "w") do file
        println(file, "X = $X")
        println(file, "Y = $Y")
    end
end

# transform_csv("divorce.csv")
# transform_csv("higher_education.csv")
transform_csv("accent.csv", true)
