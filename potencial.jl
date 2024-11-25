using CSV, DataFrames, Loess, Plots, LsqFit, Statistics, Distributions, LaTeXStrings

# Tipografía para coincidir las gráficas con el documento
Plots.default(fontfamily = ("Computer Modern"))

# Cargar datos desde el archivo CSV
df = CSV.read("potencialprueba.csv", DataFrame)


function modelchild(x, a)
	return a[1] .* x
end

# Calcular la pendiente y el intercepto
function line_through_points(x1, y1, x2, y2)
	m = (y2 - y1) / (x2 - x1)  # Pendiente
	b = y1 - m * x1            # Intercepto
	return m, b
end

# Cálculo de la primera derivada usando diferencias centrales de cuarto orden
function central_diff_first_derivative(x, y)
	n = length(x)
	h = x[2] - x[1]  # Supone espaciado uniforme
	df = zeros(n)

	for i in 3:n-2
		df[i] = (-y[i+2] + 8 * y[i+1] - 8 * y[i-1] + y[i-2]) / (12 * h)
	end

	return df[3:end-2]
end

# Cálculo de la segunda derivada usando diferencias centrales de cuarto orden
function central_diff_second_derivative(x, y)
	n = length(x)
	h = x[2] - x[1]  # Supone espaciado uniforme
	d2f = zeros(n)

	for i in 3:n-2
		d2f[i] = (-y[i+2] + 16 * y[i+1] - 30y[i] + 16y[i-1] - y[i-2]) / (12 * h^2)
	end
	return d2f[3:end-2]
end

# Cálculo de la tercera derivada usando diferencias centrales de cuarto orden
function central_finite_diff_third_derivative_4th_order(x, y)
	n = length(x)
	h = x[2] - x[1]  # Supone espaciado uniforme
	d3f = zeros(n)

	for i in 4:n-3
		d3f[i] = (-y[i+3] + 8 * y[i+2] - 13 * y[i+1] + 13 * y[i-1] - 8 * y[i-2] + y[i-3]) / (8 * h^3)
	end
	return d3f[4:end-3]
end

function find_closest_to_zero(arr::AbstractVector, tol::Float64 = 1e-6)
	# Filtrar índices donde el valor está dentro de la tolerancia
	indices = findall(x -> abs(x) <= tol, arr)

	# Retornar la primera coincidencia, si existe
	if !isempty(indices)
		idx = indices[1]  # Tomar el primer índice encontrado
		return (idx, arr[idx])  # Retornar índice y valor
	else
		return nothing  # No se encontró ningún valor dentro de la tolerancia
	end
end

# Función para calcular el error de lectura de voltaje (20V DC)
function error_lectura_voltaje(lectura::Float64, resolucion::Float64 = 0.01)
	# Error de lectura para 20V en escala DC, 1% de la lectura + 2 dígitos (resolución)
	error = lectura * 0.008 + 3 * resolucion
	return error
end

# Función para calcular el error de lectura de corriente (2000mA DC)
function error_lectura_corriente(lectura::Float64, porcentaje::Float64 = 0.001, resolucion::Float64 = 0.01)
	# Error de lectura para 2000mA en escala DC, 1% de la lectura + 2 dígitos (resolución)
	error = lectura * porcentaje + 2 * resolucion
	return error
end

function incertidumbrevel(lectura, error)

	return (3 / 2) * lectura^(1 / 2) * error

end

# Función para la línea recta
line(x) = m * x + b

function uncertantycalculator(array)

	resultado = []  # Crear una copia para guardar los resultados

	for ele in array
		# Contar la cantidad de decimales
		partes = split(string(ele), ".")

		if length(partes) == 2 && length(partes[2]) == 2
			err = error_lectura_corriente(ele, 0.001, 0.01)
			push!(resultado, err)

		elseif length(partes) == 2 && length(partes[2]) == 1
			err = error_lectura_corriente(ele, 0.002, 0.1)
			push!(resultado, err)

		elseif (ele <= 20.0)
			err = error_lectura_corriente(ele, 0.001, 0.01)
			push!(resultado, err)

		elseif (ele > 20.0)
			err = error_lectura_corriente(ele, 0.002, 0.1)
			push!(resultado, err)
		end
	end
	return resultado
end

function find_first_index(array, value)
    for (i, ele) in pairs(array)
        if ele == value
            return i  # Devuelve el índice si hay coincidencia
        end
    end
    return nothing  # Devuelve `nothing` si no se encontró el valor
end


# Extraer las columnas x e y como vectores
x = filter(!ismissing, df.V_10)
y = filter(!ismissing, df.I_10)


unx = incertidumbrevel.(x, error_lectura_voltaje.(x))

uny = uncertantycalculator(y) ./ 1000

errx = [ex * ones(2) for ex in unx]
erry = [ey * ones(2) for ey in uny]

x = x .^ (3 / 2)
y = y ./ 1000

a0 = [0.0]
# Se realiza el ajuste utilizando el algoritmo de Levenberg-Marquardt
fit_result = curve_fit(modelchild, x[1:15], y[1:15], a0)

a_fit = fit_result.param[1]
println("Coeficiente ajustado a: ", a_fit)

y_fit = modelchild(x[1:15], fit_result.param)

# Calcular el coeficiente de determinación R^2
ss_res = sum((y[1:15] .- y_fit) .^ 2)
ss_tot = sum((y[1:15] .- mean(y)) .^ 2)
r_squared = 1 - (ss_res / ss_tot)
println("Coeficiente de determinación R^2 : ", r_squared)

# Dos puntos dados
x1, y1 = x[end-7], y[end-7]
x2, y2 = x[end], y[end]


# Obtener la ecuación de la línea
m, b = line_through_points(x1, y1, x2, y2)

println("Ecuación de la línea: y = $m*x + $b")

# Ajustar el modelo LOESS

# Generar una colección equidistante de puntos
grid = range(minimum(x), stop = maximum(x), length = 100)  # 100 puntos equidistantes

grid2 = range(36, 42, length = 100)

datachild = modelchild(grid, fit_result.param)

model = loess(x, y, span = 0.58)

# Evaluar el modelo en los puntos equidistantes
smoothed_grid = predict(model, collect(grid))  # Usar predict en lugar de llamar al modelo directamente

# Graficar
img = plot(x, y,
	label = "Medición 9",
	color = :lightblue,
	alpha = 0.75,
	xerror = hcat(errx...),
	yerror = hcat(erry...),
	xlabel = L"$Voltaje^{\frac{3}{2}}\ \ [V^{\frac{3}{2}}]$", ylabel = "Corriente [A]",
	title = "Corriente a través del anado con respecto a V linealizado por ley L-C.",
)

plot!(img, grid, [smoothed_grid datachild], label = ["Curva LOESS evaluada" L"Ajuste zona L-C, $R^{2}\ =\ %$(round(r_squared, digits = 4)) $"], color = [:red :orange], alpha = [0.6 1], lw = 1.5)
plot!(img, grid2, line.(grid2), label = "Región de Saturación.", color = :pink, lw = 1.5, dpi = 620)

display(img)


newcurvedatax = x[14:end-7]
newcurvedatay = y[14:end-7]

newgrid = minimum(newcurvedatax):0.001:maximum(newcurvedatax)  # 100 puntos equidistantes 

kneemodel = loess(newcurvedatax, newcurvedatay, span = 1.0)

# Evaluar el modelo en los puntos equidistantes
kneesmoothed_grid = predict(kneemodel, collect(newgrid))  # Usar predict en lugar de llamar al modelo directamente

plot!(img, newgrid, kneesmoothed_grid, label = "Ajuste rodilla LOESS", color = :green, lw = 2, dpi = 620)

display(img)

d1f = central_diff_first_derivative(newgrid, kneesmoothed_grid)
d2f = central_diff_second_derivative(newgrid, kneesmoothed_grid)
d3f = central_finite_diff_third_derivative_4th_order(newgrid, kneesmoothed_grid)


df2domain = newgrid[3:end-2]
df3domain = newgrid[4:end-3]

pares = find_closest_to_zero(d2f, 1e-5)
display(pares[1][1])

Vp32 = df2domain[pares[1][1]]
display(Vp32)
display(df3domain[pares[1][1]-1])


vline!(img, [Vp32],
	label = L"${V_{p}}^{\frac{3}{2}}\ =\ %$(round(Vp32, digits = 3))\ V^{\frac{3}{2}} $",
	color = :gray,
	linestyle = :dash,
	legend = :outerbottom,        # Coloca la leyenda debajo
	legendcolumns = 3,
	size = (800, 600),
	margin = 5Plots.mm,          # Organiza la leyenda en 3 columnas
)

display(img)


# Graficar las derivadas
img2 = plot(newgrid[3:end-2], d1f, label = "f'",
	color = :blue,
	lw = 1.5,
	title = " Primer y segunda derivada numerica, medición 9.")
plot!(img2, newgrid[3:end-2], d2f, label = "f''", color = :purple, lw = 1.5)

vline!(img2, [Vp32],
	label = L"${V_{p}}^{\frac{3}{2}}\ =\ %$(round(Vp32, digits = 3))\ V^{\frac{3}{2}} $",
	color = :gray,
	linestyle = :dash,
	legend = :outerbottom,        # Coloca la leyenda debajo
	legendcolumns = 3,
	size = (800, 600),
	margin = 5Plots.mm, 
	dpi=640         # Organiza la leyenda en 3 columnas
)


display(img2)


img3 = plot(newgrid[4:end-3], d3f, label = "f'''",
	color = :green,
	lw = 1.5,
	title = " Tercer derivada numerica, medición 9.")

vline!(img3, [Vp32],
	label = L"${V_{p}}^{\frac{3}{2}}\ =\ %$(round(Vp32, digits = 3))\ V^{\frac{3}{2}} $",
	color = :gray,
	linestyle = :dash,
	legend = :outerbottom,        # Coloca la leyenda debajo
	size = (800, 600),
	margin = 5Plots.mm,  
	dpi=640         # Organiza la leyenda en 3 columnas
)
display(img3)

#savefig(img, "pot9.png")
#savefig(img2, "potderps9.png")
#savefig(img3, "potderst9.png")
