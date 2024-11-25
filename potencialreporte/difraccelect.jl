using ForwardDiff
using DataFrames  # Para construir dataframes
using CSV# Para escribir archivos CSV
using Plots
using LaTeXStrings
using Statistics

# Definimos las constantes
h = 6.62607015e-34  # constante de Planck
m = 9.10938356e-31  # masa del electrón
q = 1.602176634e-19  # carga del electrón

d1 = 1.23e-8 / 100
d2 = 2.13e-8 / 100

L = 0.26
deltaV = ones(12) .* 100

n = 1

# Tipografía para coincidir las gráficas con el documento
Plots.default(fontfamily=("Computer Modern"))

# Se lee el archivo CSV y se almacena en un DataFrame
data = CSV.read("datosrads.csv", DataFrame)

V = data.Voltaje .* 1000
Rextprom = [data[1:12, "Rextprom$(i)"] ./ 100 for i in 1:3]
Rintprom = [data[1:12, "Rintprom$(i)"] ./ 100 for i in 1:3]
deltaRextprom = [data[1:12, "deltaRextprom$(i)"] ./ 100 for i in 1:3]
deltaRintprom = [data[1:12, "deltaRintprom$(i)"] ./ 100 for i in 1:3]

# Funciones que calculan λ_star y d_star
function λ_star(d, R)
    return (d * R) / (n * L)
end

function λ_star_err(d, dR)
    return sqrt(((dR * d) / (n * L))^2)
end

function d_star(R, V)
    return (n * L * h) / (R * sqrt(2 * m * q * V))
end

function d_star_err(dR, R, dV, V)
    return sqrt(((-dR * h * L * n) / (sqrt(2) * R^2 * sqrt(m * q * V)))^2 + ((-dV * h * L * m * n * q) / (2 * sqrt(2) * R * (m * q * V)^(3 / 2)))^2)
end



function λ_star2(d, a)
    return (2 * d / n) * sin(a / 2)
end

function λ_star2_err(d, a, da)
    return sqrt.(((da .* d .* cos.(a / 2)) ./ n) .^ 2)  # Use . to apply element-wise
end


function d_star2(a, V)
    return (n * h) / (sqrt(2 * m * q * V) * (2 * sin(a / 2)))
end

function d_star2_err(da, a, dV, V)
    return sqrt(((-da * n * h * cot(a / 2) * csc(a / 2)) / (4 * sqrt(2) * sqrt(m * q * V)))^2 + ((-dV * h * m * n * q * csc(a / 2)) / (4 * sqrt(2) * (m * q * V)^(3 / 2)))^2)
end


function λ_boglie(V)
    return h / sqrt(2 * m * q * V)
end

function λ_boglie_err(dV, V)
    return sqrt(((-dV * h * m * q) / (2 * sqrt(2) * (m * q * V)^(3 / 2)))^2)
end


#Crear un DataFrame para almacenar los resultados
results = DataFrame()

results[!, "V"] = V

λboglie = λ_boglie.(V)
λboglie_err = λ_boglie_err.(deltaV, V)

results[!, "λ_boglie"] = λboglie
results[!, "λ_boglie_err)"] = λboglie_err

# Iterar sobre las series de resistencias y calcular los resultados
for i in 1:3
    Rext = Rextprom[i]
    Rint = Rintprom[i]
    dRext = deltaRextprom[i]
    dRint = deltaRintprom[i]

    # Calcular λ_star y su error para Rext e Rint
    λ_ext = λ_star.(d1, Rext)
    λ_int = λ_star.(d2, Rint)
    λ_ext_err = λ_star_err.(d2, dRext)
    λ_int_err = λ_star_err.(d1, dRint)

    # Calcular d_star y su error para Rext e Rint
    d_ext = d_star.(Rext, V)
    d_int = d_star.(Rint, V)
    d_ext_err = d_star_err.(dRext, Rext, deltaV, V)
    d_int_err = d_star_err.(dRint, Rint, deltaV, V)

    # Calcular el promedio entre los elementos correspondientes
    λ_prom = ((λ_ext .+ λ_int) ./ 2) .* 1e+9

    λ_desv = desviaciones_estandar = [[std([λ_ext[i], λ_int[i]]), std([λ_ext[i], λ_int[i]])] .* 1e+9 for i in 1:length(λ_ext)]

    λb_err = [[element, element].* 1e+9 for element in λboglie_err]

    P = plot(V ./ 1000, λ_prom,
        label=L"promedio entre $\lambda_{int}, \lambda_{ext}$, serie %$(i)",
        lw=3,
        seriestype=:scatter,
        xerror=[0.1, 0.1],
        yerror=hcat(λ_desv...),
        alpha=1,
        xlabel="Voltaje [kV]",
        ylabel=" λ [nm]",
        title="λ del haz de electrones ",
        legend=:outerbottom,
        dpi=620
    )

    plot!(V ./ 1000, λboglie .* 1e+9,
        label="λ hip. de Broglie",
        lw=3,
        xerror=[0.1, 0.1],
        yerror=hcat(λb_err...),
        seriestype=:scatter,
        dpi=620
    )

    savefig(P, "lambdapm$(i).png")
    display(P)

    d_ex_e = [[element, element] for element in d_ext_err]
    d_in_e = [[element, element] for element in d_int_err]

    D = plot(V ./ 1000, [(d_ext.*1e+9) (d_int.*1e+9) ],
        label=[L"$d^{*}$ anillo ext serie %$(i)" L"$d^{*}$ anillo int serie %$(i)"],
        lw=[3 3],
        seriestype=[:scatter :scatter],
        xerror=[0.1, 0.1],
        yerror=[hcat(d_ex_e...) hcat(d_in_e...)],
        alpha=1,
        xlabel="Voltaje [kV]",
        ylabel=L" $d^{*}$ [nm]",
        title=L"$d^{*}$ para el anillo interno y externo",
        legend=:outerbottom,
        dpi=620
    )

    savefig(D, "dpm$(i).png")
    display(D)

    errdextd1 = [ (abs(element - d1)/d1)*100 for element in d_ext] 
    errdintd2 = [ (abs(element - d2)/d2)*100 for element in d_int] 

    E = plot(V ./ 1000, [errdintd2  errdextd1 ],
        label=[L" $d^{*}$ anillo interno serie %$(i), $d_{2}$" L"$d^{*}$ anillo externo serie %$(i),  $d_{1}$"],
        lw=[3 3],
        seriestype=[:scatter :scatter],
        alpha=1,
        xlabel="Voltaje [kV]",
        ylabel="error %",
        title=L" Error entre $d^{*}$, d para el anillo interno y externo",
        legend=:outerbottom,
        dpi=620
    )

    savefig(E, "errdpm$(i).png")
    display(E)
    
    # Agregar columnas al DataFrame de resultados
    results[!, "λ_ext_series$(i)"] = λ_ext
    results[!, "λ_ext_err_series$(i)"] = λ_ext_err
    results[!, "λ_int_series$(i)"] = λ_int
    results[!, "λ_int_err_series$(i)"] = λ_int_err
    results[!, "d_ext_series$(i)"] = d_ext
    results[!, "d_ext_err_series$(i)"] = d_ext_err
    results[!, "d_int_series$(i)"] = d_int
    results[!, "d_int_err_series$(i)"] = d_int_err
    
end



# Guardar los resultados en un archivo CSV
CSV.write("resultados_con_errores_primermet.csv", results)

angs = CSV.read("angsvrad.csv", DataFrame)

V = angs.V
deltaV = ones(6) * 100
dosthetaint = angs.aprim
dosthetainterr = angs.deltaaprim
dosthetaext = angs.asec
dosthetaexterr = angs.deltaasec

angsresults = DataFrame()

angsresults[!, "V"] = V

int_d_star2_vals = d_star2.(dosthetaint, V)
int_d_star2_errs = d_star2_err.(dosthetainterr, dosthetaint, deltaV, V)

ext_d_star2_vals = d_star2.(dosthetaext, V)
ext_d_star2_errs = d_star2_err.(dosthetaexterr, dosthetaext, deltaV, V)

int_λ_star2_vals = λ_star2.(d2, dosthetaint)
int_λ_star2_errs = λ_star2_err(d2, dosthetaint, dosthetainterr)

ext_λ_star2_vals = λ_star2.(d1, dosthetaext)
ext_λ_star2_errs = λ_star2_err(d1, dosthetaext, dosthetaexterr)

λ_prom = ((int_λ_star2_vals .+ ext_λ_star2_vals) ./ 2) .* 1e+9

λ_desvm = [std([int_λ_star2_vals[i],ext_λ_star2_vals[i]])* 1e+9 for i in 1:length(int_λ_star2_vals)]
λb_err = λboglie_err[1:6].*1e+9

W = plot(V ./ 1000, [λ_prom (λboglie[1:6] .* 1e+9)],
        label= [L"promedio entre $\lambda_{int}, \lambda_{ext}$, segundo met" "λ hip. de Broglie"],
        lw=[3 3],
        seriestype=[:scatter :scatter],
        xerror=[0.1 0.1],
        yerror=[λ_desvm λb_err],
        alpha=1,
        xlabel="Voltaje [kV]",
        ylabel=" λ [nm]",
        title="λ del haz de electrones ",
        legend=:outerbottom,
        dpi=620
    )

    savefig(W, "lambdasm.png")
    display(W)

    d_ex_e = [[element, element] for element in ext_d_star2_errs]
    d_in_e = [[element, element] for element in int_d_star2_errs]

    D = plot(V ./ 1000, [(ext_d_star2_vals.*1e+9) (int_d_star2_vals.*1e+9) ],
        label=[L"$d^{*}$ anillo ext segundo met" L"$d^{*}$ anillo int segundo met"],
        lw=[3 3],
        seriestype=[:scatter :scatter],
        xerror=[0.1, 0.1],
        yerror=[hcat(d_ex_e...) hcat(d_in_e...)],
        alpha=1,
        xlabel="Voltaje [kV]",
        ylabel=L" $d^{*}$ [nm]",
        title=L"$d^{*}$ para el anillo interno y externo",
        legend=:outerbottom,
        dpi=620
    )

    savefig(D, "dsm.png")
    display(D)

    errdextd1 = [ (abs(element - d1)/d1)*100 for element in ext_d_star2_vals] 
    errdintd2 = [ (abs(element - d2)/d2)*100 for element in int_d_star2_vals] 

    E = plot(V ./ 1000, [errdintd2  errdextd1 ],
        label=[L" $d^{*}$ anillo interno segundo met, $d_{2}$" L"$d^{*}$ anillo externo segundo met,  $d_{1}$"],
        lw=[3 3],
        seriestype=[:scatter :scatter],
        alpha=1,
        xlabel="Voltaje [kV]",
        ylabel="error %",
        title=L" Error entre $d^{*}$, d para el anillo interno y externo",
        legend=:outerbottom,
        dpi=620
    )

    savefig(E, "errdsm.png")
    display(E)


angsresults[!, "d* int"] = int_d_star2_vals
angsresults[!, "d* int err"] = int_d_star2_errs
angsresults[!, "d* ext"] = ext_d_star2_vals
angsresults[!, "d* ext err"] = ext_d_star2_errs
angsresults[!, "λ* int"] = int_λ_star2_vals
angsresults[!, "λ* int err"] = int_λ_star2_errs
angsresults[!, "λ* ext"] = ext_λ_star2_vals
angsresults[!, "λ* ext err"] = ext_λ_star2_errs

#CSV.write("resultados_con_errores_secmet.csv", angsresults)