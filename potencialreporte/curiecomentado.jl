#=--------------------------------------
# Perez Julio
# Azpeitia Alejandro
# Sanchez Fabricio

# Análisis Curie
----------------------------------------=#

# Importar librerías necesarias
using CSV
using DataFrames
using Plots
using LaTeXStrings
using Statistics
using LinearAlgebra

# Configuración de la tipografía para que coincida con el documento
Plots.default(fontfamily=("Computer Modern"))

# Leer los archivos CSV con los datos experimentales
df1 = CSV.read("curie/curierfriocalienteuno.csv", DataFrame)
df2 = CSV.read("curie/curiercalientefriouno.csv", DataFrame)
df3 = CSV.read("curie/curierfriocalientedos.csv", DataFrame)
df4 = CSV.read("curie/curiercalientefriodos.csv", DataFrame)

#=
# Función para combinar filas duplicadas en un DataFrame basándose en una columna de interés
# Agrupa las filas por un valor en la columna y calcula el promedio de las demás columnas
function merge_duplicate_rows!(df::DataFrame, col_name::Symbol)
        # Obtener los valores únicos en la columna de interés
        unique_values = unique(df[!, col_name])
        # Inicializar un DataFrame vacío con las mismas columnas que `df`
        clean_df = DataFrame(df[1:0, :])

        for element in unique_values
            # Encontrar los índices de las filas donde el valor coincide con el elemento actual
            matching_indices = findall(isequal(element), df[!, col_name])

            if length(matching_indices) > 0
                # Calcular el promedio para cada columna para filas con el mismo valor
                averaged_row = [mean(df[i, j] for i in matching_indices) for j in 1:ncol(df)]

                # Agregar la fila promediada al DataFrame limpio
                push!(clean_df, averaged_row, promote=true)
            end
        end

        return clean_df
    end

# Aplicar la función para combinar filas duplicadas en cada DataFrame
df1 = merge_duplicate_rows!(df1,:Temp)
df2 = merge_duplicate_rows!(df2,:Temp)
df3 = merge_duplicate_rows!(df3,:Temp)
df4 = merge_duplicate_rows!(df4,:Temp)

# Guardar los DataFrames combinados en nuevos archivos CSV
CSV.write("curierfriocalienteuno.csv", df1)
CSV.write("curiercalientefriouno.csv", df2)
CSV.write("curierfriocalientedos.csv", df3)
CSV.write("curiercalientefriodos.csv", df4)
=#

# Extraer los datos de interés de voltaje (VS) y temperatura (T)
VS = [Float64.(df1.VRMSsec), Float64.(df2.VRMSsec), Float64.(df3.VRMSsec), Float64.(df4.VRMSsec)]
T = [Float64.(df1.Temp), Float64.(df2.Temp), Float64.(df3.Temp), Float64.(df4.Temp)]
VSA = [Float64.(df1.VRMSnoamosec), Float64.(df2.VRMSnoamosec), Float64.(df3.VRMSnoamosec), Float64.(df4.VRMSnoamosec)]

# Función para suavizar datos usando un ajuste polinómico local (LOESS)
function loess(x, y, h; degree=2)
    n = length(x)  # Número de puntos de datos
    smoothed_y = zeros(n)  # Vector para almacenar los valores suavizados

    for i in 1:n
        # Definir el vecindario local
        d = abs.(x .- x[i])
        weights = exp.(-d .^ 2 / (h^2))  # Función de peso gaussiana

        # Preparar matrices locales para el ajuste
        X = hcat(ones(n), x, x .^ 2)  # Matriz de diseño para polinomios cuadráticos
        W = Diagonal(weights)

        # Ajuste polinomial local utilizando mínimos cuadrados ponderados
        A = X' * W * X
        b = X' * W * y

        if det(A) == 0
            # Manejar caso de matriz singular
            coefs = [0.0, 0.0, 0.0]
        else
            coefs = A \ b
        end

        # Evaluar el polinomio en x[i]
        smoothed_y[i] = coefs[1] + coefs[2] * x[i] + coefs[3] * x[i]^2
    end

    return smoothed_y
end

# Función para realizar interpolación de Lagrange
function lagrange_interpolation(x_points, y_points, x)
    n = length(x_points)
    P_x = 0.0

    for i in 1:n
        L_i = 1.0
        for j in 1:n
            if i != j
                L_i *= (x - x_points[j]) / (x_points[i] - x_points[j])
            end
        end
        P_x += y_points[i] * L_i
    end

    return P_x
end

# Cálculo de la primera derivada usando diferencias centrales de cuarto orden
function central_diff_first_derivative(x, y)
    n = length(x)
    h = x[2] - x[1]  # Supone espaciado uniforme
    df = zeros(n)

    for i in 3:n-2
        df[i] = (-y[i+2] + 8 * y[i+1] - 8 * y[i-1] + y[i-2]) / (12 * h)
    end

    return df
end

# Cálculo de la segunda derivada usando diferencias centrales de cuarto orden
function central_diff_second_derivative(x, y)
    n = length(x)
    h = x[2] - x[1]  # Supone espaciado uniforme
    d2f = zeros(n)

    for i in 3:n-2
        d2f[i] = (-y[i+2] + 16 * y[i+1] - 30y[i] + 16y[i-1] - y[i-2]) / (12 * h^2)
    end
    return d2f
end

# Cálculo de la tercera derivada usando diferencias centrales de cuarto orden
function central_finite_diff_third_derivative_4th_order(x, y)
    n = length(x)
    h = x[2] - x[1]  # Supone espaciado uniforme
    d3f = zeros(n)

    for i in 4:n-3
        d3f[i] = (-y[i+3] + 8 * y[i+2] - 13 * y[i+1] + 13 * y[i-1] - 8 * y[i-2] + y[i-3]) / (8 * h^3)
    end
    return d3f
end

# Parámetro de suavizado
h = 5
domain = 6:0.25:33  # Dominio para interpolación

# Suavizar los datos de voltaje usando LOESS
smoothed_Vs = [loess(T[i], VS[i], h) for i in 1:4]

# Realizar interpolación de Lagrange
interpolated_Vs = [[lagrange_interpolation(T[i], smoothed_Vs[i], x) for x in domain] for i in 1:4]

# Cálculo de derivadas
d1f = [central_diff_first_derivative(domain, interpolated_Vs[i]) for i in 1:4]
d2f = [central_diff_second_derivative(domain, interpolated_Vs[i]) for i in 1:4]
d3f = [central_finite_diff_third_derivative_4th_order(domain, interpolated_Vs[i]) for i in 1:4]

# Lista para almacenar el punto de inflexión mínimo para cada ciclo
dfmin = []

# Encontrar el mínimo de la primera derivada para cada ciclo
for i in 1:4
    dfmini = first([(j, t - 0.5, d1f[1][j]) for (j, t) in enumerate(domain[3:length(domain)-2]) if d1f[i][j] == minimum(d1f[i])])
    push!( dfmin, dfmini)
end

Ts = [ dfmin[i][2] for i in 1:4]

display(mean(Ts))
display(std(Ts))

y_err =  [(v * 0.012 + 0.005) * ones(2) for v in VS[1]]

P = plot(T[1], VS[1],
    label="Datos Medidos",
    lw=3,
    seriestype=:scatter,
    xerror=[0.25, 0.25],
    yerror=hcat(y_err...),
    alpha=1,
    xlabel="Temperatura [C°]",
    ylabel="Voltaje [V]",
    title=L"$V_{rms}$ secundario, ciclo n = 1",
    dpi=620
)

plot!(P, T[1], smoothed_Vs[1],
    label="loess",
    lw=3,
    seriestype=:line,
    alpha=0.6,
    dpi=620
)

plot!(P, domain, interpolated_Vs[1],
    label="interpolación Lagrande",
    lw=1.5, alpha=0.75)

vline!([dfmin[1][2]], label="punto de inflexión", linestyle=:dash, alpha=0.5, color=:gray)
display(P)
savefig(P, "VRMSsec1.png")


for i in 2:4


    yerr =  [(sqrt( ((8e-5*1e+5)/8.4e+4)^2 + ((1000*VSA[i][j])/8.4e+4)^2 + ((776*1e+3*VSA[i][j])/(8.4e+4)^2)^2)) * ones(2) for j in 1:length(VS[i]) ]
   

    P = plot(T[i], VS[i],
        label="Datos Medidos",
        lw=3,
        seriestype=:scatter,
        xerror=[0.25, 0.25],
        yerror=hcat(y_err...),
        alpha=1,
        xlabel="Temperatura [C°]",
        ylabel="Voltaje [V]",
        title=L"$V_{rms}$ secundario, ciclo n = %$(i)",
        dpi=620
    )

    plot!(P, T[i], smoothed_Vs[i],
        label="loess",
        lw=3,
        seriestype=:line,
        alpha=0.6,
        dpi=620
    )

    plot!(P, domain, interpolated_Vs[i],
        label="interpolación Lagrande",
        lw=1.5, alpha=0.75)

    vline!([dfmin[i][2]], label="punto de inflexión", linestyle=:dash, alpha=0.5, color=:gray)
    display(P)
    savefig(P, "VRMSsec$(i).png")

end

#gráficar las fer
for i in 1:4

    d1 = d1f[i]
    d2 = d2f[i]
    d3 = d3f[i]

    P = plot(domain, [d1 d2 d3],
        title="Derivadas, ciclo n = $(i)",
        xlabel="Temperatura [C°]",
        label=[L"$V^{\prime}(T)\ [VC^{\circ^{-1}}]$" L"$V^{\prime \prime}(T)\ [VC^{\circ^{-2}}]$" L"$V^{\prime \prime \prime}(T)\ [VC^{\circ^{-3}}]$"],
        dpi=620
    )

    vline!([dfmin[i][2]], label="punto de inflexión", linestyle=:dash, alpha=0.5, color=:gray)
    hline!([0], label=false, linestyle=:dash, alpha=0.5, color=:gray)

    display(P)
    savefig(P, "derivadas$(i).png")

end



