#=Perez Flores Julio A.
  Analisis de datos practica Faraday Lab. Contempo I.
=#

using CSV # para cargar archivos csv
using Plots # para gráficar
using DataFrames # para convertir un CSV en DataFrame
using LsqFit # Ajustes no lineales
using Statistics #usar promedio
using LaTeXStrings # Para escribir mates
using LinearAlgebra # manipular matrices

# Tipografía para coincidir las gráficas con el documento
Plots.default(fontfamily=("Computer Modern"))

# Se lee el archivo CSV y se almacena en un DataFrame
data = CSV.read("faradaySI.csv", DataFrame)


#= 
------------------------------------------------
ESTIMACIÓN MAGNITUD CAMPO B con datasheet
------------------------------------------------
=#

TDS = [0*0.0001 ,3000*0.0001, 6000*0.0001, 7000*0.0001]
ADS = [0 , 2, 4, 5]

# Se define la función de modelo (lineal)
function model(x, a)
    return a[1] .* x
end

a = [0.0]
Bfit = curve_fit(model, ADS, TDS, a)

AE = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,  5.0]
BE = model(AE, Bfit.param) #B's para las A's usadas

#= 
------------------------------------------------
PARTE EXPERIMENTAL
------------------------------------------------
=#

#Magnitud Campo B registrada SF-2
BN = [data[1:11, "BN1"],
    data[1:11, "BN2"],
    data[1:11, "BN3"]]

#Error del campo B medido * l del vidrio
errBIM1 = []

for i in 1:3
    er= []
    for j in 2:11
        dl = 0.0505*0.00005 + BN[i][j]*0.0001
        push!(er, [dl, dl])
    end
    push!(errBIM1, er)
end

#theta registrada SF-2
theta = [data[1:11, "°N1"] .* (pi / 180),
    data[1:11, "°N2"] .* (pi / 180),
    data[1:11, "°N3"] .* (pi / 180)]


for i in 1:3
    local a0 = [0.0]
    # Se realiza el ajuste utilizando el algoritmo de Levenberg-Marquardt
    fit_result = curve_fit(model, BN[i] .* 0.0505, theta[i], a0)

    # Se calcula la función de ajuste para los datos
    y_fit = model(BN[i] .* 0.0505, fit_result.param)

    # Se extrae la matriz de covarianza
    cov_matrix = estimate_covar(fit_result)

    # Se extraen los elementos diagonales (varianzas) de la matriz de covarianza
    variances = diag(cov_matrix)

    # Se calculan los errores estándar como la raíz cuadrada de las varianzas
    std_errors = sqrt.(variances)

    # Se imprimen los errores estándar de los coeficientes
    for (j, std_err) in enumerate(std_errors)
        println("Error estándar del coeficiente muestra $(i), coeficiente $(j): ", std_err)
    end

    # Se imprimen los valores de los coeficientes del ajuste
    a_fit = fit_result.param[1]
    println("Coeficiente ajustado a (muestra $(i)): ", a_fit)

    # Calcular el coeficiente de determinación R^2
    ss_res = sum((theta[i] .- y_fit) .^ 2)
    ss_tot = sum((theta[i] .- mean(theta[i])) .^ 2)
    r_squared = round(1 - (ss_res / ss_tot), digits=3)
    println("Coeficiente de determinación R^2 (muestra $(i)): ", r_squared)

    # Crear la gráfica
    img = plot(
        BN[i] .* 0.0505, theta[i],
        yerror=0.00872665,
        xerror=hcat(errBIM1[i]...),
        label="SF-02, l = 0.0505 m, medición: $(i)",
        seriestype=:scatter,
        ylabel=L"$\Delta \theta \ [rad]$",
        xlabel=L"$ B*l \ [Tm] $",
        title=L"Relación $Bl,\ \Delta \theta$",
        linecolor=:blue,
        dpi=620
    )

    # Se agrega la curva de ajuste lineal
    x = 0.0001*0.0505:0.00001:0.04*0.0505
    plot!(img, x, model(x, fit_result.param),
        linecolor=:purple,
        label=L"ajuste lineal $R^2\ =\ %$(r_squared)$ 
                $V^{*}  = %$(round(a_fit, digits =2)),\ \delta V^{*}  = %$(round(std_errors[1], digits =3)) $"
    )

    # Mostrar la gráfica
    display(img)

    # Exportar la gráfica
    savefig(img, "faradaymat1med$(i).png")
end

#Magnitud Campo B registrada SF-5
BNM2 = [data[1:11, "BN1M2"],
    data[1:11, "BN2M2"],
    data[1:11, "BN3M2"]]

#error de B*l para SF-5
errBIM2 = []

#Error del campo B medido * l del vidrio
for i in 1:3
    er= []
    for j in 2:11
        dl = 0.0404*0.00005 + BNM2[i][j]*0.0001
        push!(er, [dl, dl])
    end
    push!(errBIM2, er)
end

#theta registrada SF-5
thetaM2 = [data[1:11, "°N1M2"] .* (pi / 180),
    data[1:11, "°N2M2"] .* (pi / 180),
    data[1:11, "°N3M2"] .* (pi / 180)]


for i in 1:3
    local a0 = [0.0]

    # Se realiza el ajuste utilizando el algoritmo de Levenberg-Marquardt
    fit_result = curve_fit(model, BNM2[i] .* 0.0404, thetaM2[i], a0)

    # Se calcula la función de ajuste para los datos
    y_fit = model(BNM2[i] .* 0.0404, fit_result.param)

    # Se extrae la matriz de covarianza
    cov_matrix = estimate_covar(fit_result)

    # Se extraen los elementos diagonales (varianzas) de la matriz de covarianza
    variances = diag(cov_matrix)

    # Se calculan los errores estándar como la raíz cuadrada de las varianzas
    std_errors = sqrt.(variances)

    # Se imprimen los errores estándar de los coeficientes
    for (j, std_err) in enumerate(std_errors)
        println("Error estándar del coeficiente muestra $(i), coeficiente $(j): ", std_err)
    end

    # Se imprimen los valores de los coeficientes del ajuste
    a_fit = fit_result.param[1]
    println("Coeficiente ajustado a (muestra $(i)): ", a_fit)

    # Calcular el coeficiente de determinación R^2
    ss_res = sum((thetaM2[i] .- y_fit) .^ 2)
    ss_tot = sum((thetaM2[i] .- mean(thetaM2[i])) .^ 2)
    r_squared = round(1 - (ss_res / ss_tot), digits=3)
    println("Coeficiente de determinación R^2 (muestra $(i)): ", r_squared)

    # Crear la gráfica
    img = plot(
        BNM2[i] .* 0.0404, thetaM2[i],
        yerror=0.00872665,
        xerror=hcat(errBIM2[i]...),
        label="SF-05, l = 0.0404 m, medición: $(i)",
        seriestype=:scatter,
        ylabel=L"$\Delta \theta \ [rad]$",
        xlabel=L"$ B*l \ [Tm] $",
        title=L"Relación $Bl,\ \Delta \theta$",
        linecolor=:blue,
        dpi=620
    )

    # Se agrega la curva de ajuste lineal
    x = 0.0001*0.0404:0.00001:0.035*0.0404
    plot!(img, x, model(x, fit_result.param),
        linecolor=:purple,
        label=L"ajuste lineal $R^2\ =\ %$(r_squared)$ 
                $V^{*} = %$(round(a_fit, digits =2)),\ \delta V^{*} = %$(round(std_errors[1], digits =3)) $"
    )

    # Mostrar la gráfica
    display(img)

    # Exportar la gráfica
    savefig(img, "faradaymat2med$(i).png")
end

#=
--------------------------------------
Determinación de V* con B del datasheet
---------------------------------------
=#

# SF-2
for i in 1:3
    local a0 = [0.0]
    # Se realiza el ajuste utilizando el algoritmo de Levenberg-Marquardt
    fit_result = curve_fit(model, BE .* 0.0505, theta[i], a0)

    # Se calcula la función de ajuste para los datos
    y_fit = model(BE .* 0.0505, fit_result.param)

    # Se extrae la matriz de covarianza
    cov_matrix = estimate_covar(fit_result)

    # Se extraen los elementos diagonales (varianzas) de la matriz de covarianza
    variances = diag(cov_matrix)

    # Se calculan los errores estándar como la raíz cuadrada de las varianzas
    std_errors = sqrt.(variances)

    # Se imprimen los errores estándar de los coeficientes
    for (j, std_err) in enumerate(std_errors)
        println("Error estándar del coeficiente muestra $(i), coeficiente $(j): ", std_err)
    end

    # Se imprimen los valores de los coeficientes del ajuste
    a_fit = fit_result.param[1]
    println("Coeficiente ajustado a (muestra $(i)): ", a_fit)

    # Calcular el coeficiente de determinación R^2
    ss_res = sum((theta[i] .- y_fit) .^ 2)
    ss_tot = sum((theta[i] .- mean(theta[i])) .^ 2)
    r_squared = round(1 - (ss_res / ss_tot), digits=3)
    println("Coeficiente de determinación R^2 (muestra $(i)): ", r_squared)

    # Crear la gráfica
    img = plot(
        BE .* 0.0505, theta[i],
        yerror=0.00872665,
        label="SF-02, l = 0.0505 m, medición: $(i)",
        seriestype=:scatter,
        ylabel=L"$\Delta \theta \ [rad]$",
        xlabel=L"$ B*l \ [Tm] $",
        title=L"Relación $Bl,\ \Delta \theta$",
        linecolor=:blue,
        dpi=620
    )

    # Se agrega la curva de ajuste lineal
    x = BE[1]*0.0505:0.00001:BE[11]*0.0505
    plot!(img, x, model(x, fit_result.param),
        linecolor=:purple,
        label=L"ajuste lineal $R^2\ =\ %$(r_squared)$ 
                $V^{*}  = %$(round(a_fit, digits =2)),\ \delta V^{*}  = %$(round(std_errors[1], digits =3)) $"
    )

    # Mostrar la gráfica
    display(img)

    # Exportar la gráfica
    savefig(img, "faradaymat1medBE$(i).png")
end

#SF-5
for i in 1:3
    local a0 = [0.0]

    # Se realiza el ajuste utilizando el algoritmo de Levenberg-Marquardt
    fit_result = curve_fit(model, BE .* 0.0404, thetaM2[i], a0)

    # Se calcula la función de ajuste para los datos
    y_fit = model(BE .* 0.0404, fit_result.param)

    # Se extrae la matriz de covarianza
    cov_matrix = estimate_covar(fit_result)

    # Se extraen los elementos diagonales (varianzas) de la matriz de covarianza
    variances = diag(cov_matrix)

    # Se calculan los errores estándar como la raíz cuadrada de las varianzas
    std_errors = sqrt.(variances)

    # Se imprimen los errores estándar de los coeficientes
    for (j, std_err) in enumerate(std_errors)
        println("Error estándar del coeficiente muestra $(i), coeficiente $(j): ", std_err)
    end

    # Se imprimen los valores de los coeficientes del ajuste
    a_fit = fit_result.param[1]
    println("Coeficiente ajustado a (muestra $(i)): ", a_fit)

    # Calcular el coeficiente de determinación R^2
    ss_res = sum((thetaM2[i] .- y_fit) .^ 2)
    ss_tot = sum((thetaM2[i] .- mean(thetaM2[i])) .^ 2)
    r_squared = round(1 - (ss_res / ss_tot), digits=3)
    println("Coeficiente de determinación R^2 (muestra $(i)): ", r_squared)

    # Crear la gráfica
    img = plot(
        BE .* 0.0404, thetaM2[i],
        yerror=0.00872665,
        label="SF-05, l = 0.0404 m, medición: $(i)",
        seriestype=:scatter,
        ylabel=L"$\Delta \theta \ [rad]$",
        xlabel=L"$ B*l \ [Tm] $",
        title=L"Relación $Bl,\ \Delta \theta$",
        linecolor=:blue,
        dpi=620
    )

    # Se agrega la curva de ajuste lineal
    x = BE[1]*0.0404:0.00001:BE[11]*0.0404
    plot!(img, x, model(x, fit_result.param),
        linecolor=:purple,
        label=L"ajuste lineal $R^2\ =\ %$(r_squared)$ 
                $V^{*} = %$(round(a_fit, digits =2)),\ \delta V^{*} = %$(round(std_errors[1], digits =3)) $"
    )

    # Mostrar la gráfica
    display(img)

    # Exportar la gráfica
    savefig(img, "faradaymat2medBE$(i).png")
end

#= 
------------------------------------------------
Diferencia entre I_o con B=0, a la I mas cercana a Io 
para cada cambio de B.
------------------------------------------------
=#

#Para SF-2

ITM1 = [data[1:11, "IN1"], data[1:11, "IN2"], data[1:11, "IN3"]]
errITM1 = []

#Errores para Iaprox-I_o
for i in 1:3
    er= []
    for j in 2:11
        dl = ((ITM1[i][j] * 0.05) + 2) * 2
        push!(er, [dl, dl])
    end
    push!(errITM1, er)
end


IoM1 = [ITM1[i][1] for i in 1:3]
IaM1 = [ITM1[i][2:11] for i in 1:3]

# Calcular deltaI usando broadcasting para restar cada valor de IoM1 a su correspondiente en IaM1
deltaI = [IaM1[i] .- IoM1[i] for i in 1:3]
BNO = [BN[i][2:11] for i in 1:3]

# Gráfico con los errores
deltaI_plot = plot(
    ylabel=L"$I_{0}-I_{approx}\ [Lux]$",
    xlabel=L"$ B \ [T] $",
    title=L"$I_{0}-I_{approx}$, en función de $B$",
    legend=:outerbottom,
     legendcolumns=3, 
    dpi=620
)

Ber = [[0.00005, 0.00005] for i in 1:10]

for i in 1:3
    # Gráfico con los errores
    plot!(
        deltaI_plot,
        BNO[i],
        deltaI[i],
        yerror = hcat(errITM1[i]...),
        xerror = hcat(Ber[i]...),
        ylabel=L"$I_{0}-I_{approx}\ [Lux]$",
        xlabel=L"$ B \ [T] $",
        title=L"$I_{0}-I_{approx}$, en función de $B$, para SF-2",
        label=L"n: %$(i),  $I_{o} =\ %$(IoM1[i]) \pm %$(round( (IoM1[i]*0.05)+2, digits = 0))$ Lux ",
        dpi=620
    )
end

display(deltaI_plot)
savefig(deltaI_plot, "i_iomat1.png")


#Para SF-2

ITM2 = [data[1:11, "IN1M2"], data[1:11, "IN2M2"], data[1:11, "IN3M2"]]
errITM2 = []

#Errores para Iaprox-I_o
for i in 1:3
    er= []
    for j in 2:11
        dl = ((ITM2[i][j] * 0.05) + 2) * 2
        push!(er, [dl, dl])
    end
    push!(errITM2, er)
end

IoM2 = [ITM2[i][1] for i in 1:3]
IaM2 = [ITM2[i][2:11] for i in 1:3]

# Calcular deltaI usando broadcasting para restar cada valor de IoM1 a su correspondiente en IaM1
deltaI = [IaM2[i] .- IoM2[i] for i in 1:3]
BNO = [BNM2[i][2:11] for i in 1:3]

# Gráfico con los errores
deltaI_plot = plot(
    ylabel=L"$I_{0}-I_{approx}\ [Lux]$",
    xlabel=L"$ B \ [T] $",
    title=L"$I_{0}-I_{approx}$, en función de $B$",
    legend=:outerbottom,
     legendcolumns=3, 
    dpi=620
)

#Error del campo B medido
Ber = [[0.00005, 0.00005] for i in 1:10]

for i in 1:3
    # Gráfico con los errores
    plot!(
        deltaI_plot,
        BNO[i],
        deltaI[i],
        yerror = hcat(errITM2[i]...),
        xerror = hcat(Ber[i]...),
        ylabel=L"$I_{0}-I_{approx}\ [Lux]$",
        xlabel=L"$ B \ [T] $",
        title=L"$I_{0}-I_{approx}$, en función de $B$ para SF-5",
        label=L"n: %$(i),  $I_{o} =\ %$(IoM2[i]) \pm %$(round( (IoM2[i]*0.05)+2, digits = 0))$ Lux ",
        dpi=620
    )
end

display(deltaI_plot)
savefig(deltaI_plot, "i_iomat2.png")