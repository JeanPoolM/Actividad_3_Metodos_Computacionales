using LinearAlgebra
using DataFrames
using PrettyTables
using Clustering
using Statistics
using Plots
using CSV

## CALCULAR LA MATRIZ DE ADMITANCIA NODAL

function calcular_Ynn_Yns(lines,nodes)
    """
    Entradas:   lines: DataFrame con los datos de las lineas del sistema
                nodes : DataFrame con los datos de los nodos del sistema
    Salida :    Ynn : Matriz de admitancia de los nodos diferentes al Slack.
                Yns : Matriz de admitancia de los nodos diferentes al Slack con respecto al Slack.
    """
    num_nodes = nrow(nodes)
    num_lines = nrow(lines)
    Ybus = zeros(num_nodes, num_nodes)*1im
    s = parse(Int64, string(nodes[nodes.TYPE .== 3, "NAME"][1][2]))
    for k = 1:num_lines
        # Nodo de envío
        n1 = parse(Int64,string(lines.FROM[k])[2:end])
        # Nodo de recibo
        n2 = parse(Int64,string(lines.TO[k])[2:end])
        # Admitancia de la línea
        yL = 1/(lines.R[k]+lines.X[k]*1im)
        # Susceptancia de la línea
        Bs = lines.B[k]*1im/2
        # Valor del TAP
        t = lines.TAP[k]
        if lines.TAP[k] == 0
            Ybus[n1,n1] += yL + Bs   # Dentro de la diagonal
            Ybus[n1,n2] -= yL        # Fuera de la diagonal
            Ybus[n2,n1] -= yL        # Fuera de la diagonal
            Ybus[n2,n2] += yL + Bs   # Dentro de la diagonal
        else
            Ybus[n1,n1] += (t^2 - t)*yL  # Dentro de la diagonal
            Ybus[n1,n2] -= t*yL          # Fuera de la diagonal
            Ybus[n2,n1] -= t*yL          # Fuera de la diagonal
            Ybus[n2,n2] += (1-t)*yL      # Dentro de la diagonal
        end
    end
    # Separación de la matriz Ybus en Ynn y Yns
    Ynn = Ybus[setdiff(1:end,s),setdiff(1:end,s)]
    Yns = Ybus[s+1:end,s]
    return Ynn, Yns
end

## FLUJO ESTÁTICO
function Flujo_punto_fijo(lines, nodes)
    """
    Entradas:   lines: DataFrame con los datos de las lineas del sistema
                nodes : DataFrame con los datos de los nodos del sistema
    Salida :    Vn : Vector de tensiones nodales obtenido del punto fijo.
    """
    num_nodes = nrow(nodes)
    # Cálculo de Ynn y Yns
    Ynn, Yns = calcular_Ynn_Yns(lines,nodes)
    iYnn = inv(Ynn)

    # Se inicializa Vn y Vs.
    Vn = ones(num_nodes-1) + 1im*zeros(num_nodes-1)
    Vs = 1 + 0*1im
    
  
    # Se define Sn
    Sn = (nodes.PGEN[2:end] - nodes.PLOAD[2:end]) + (nodes.QGEN[2:end]-nodes.QLOAD[2:end])*1im
    
    # Se comienza la iteración de punto fijo
    errores = zeros(4)
    ite = zeros(4)
   for i= 1:4
        # Ecuación de Vn para punto fijo.
        Vn_ant = Vn
        Vn = iYnn*(conj.(Sn./Vn) .- Yns*Vs) 
        # Se calcula el error
        error = norm(Vn - Vn_ant)
        errores[i] = error
        ite[i] = i
    end
    P = plot(ite, errores, yscale=:log10, xlabel= "Iteraciones", ylabel="Errores")
    display(P)
    # Se crea una tabla con los resultados de Vn
    col_nodos = collect(1:num_nodes-1)
    df = DataFrame(Nodos = col_nodos, Vn = Vn)
    tabla = pretty_table(df)
    display(tabla)
    return Vn    
end

lines = DataFrame(CSV.File("lines.csv"))
nodes = DataFrame(CSV.File("nodes.csv"))
Vn = Flujo_punto_fijo(lines, nodes)

## AJUSTE DE DATAFRAMES.

# Ajustar la base da datos SolarData a una matriz MuestrasxDias
SolarData_df = DataFrame(CSV.File("SolarData.csv"))
Demandas_df= DataFrame(CSV.File("Demandas.csv"))
SolarData_df2 = DataFrame()
for i = 1:288:nrow(SolarData_df)
    SolarData_df2[!, SolarData_df.Fecha[i]]= SolarData_df.Potencia[i:i+287]
end

# Se guarda el dataframe en un CSV
CSV.write("solardata_matriz.csv", SolarData_df2)

# Convertir SolarData_df2 a una matriz numérica (288 filas × 365 columnas)
matriz_solar = Matrix(SolarData_df2)

# Se gráfica las curvas de generación por día para los 365 días.
theme(:dark)
p1 = plot(size=(800, 600), legend=false)
for i in 1:365
    plot!(p1, matriz_solar[:, i], label = false, alpha=0.5)
end
title!("Potencia Generada por día para los 365 días")
xlabel!("Intervalos de 5 min")
ylabel!("Potencia Generada")

# Aplicar K-means utilizando 4 Clusters (4 Casos Típicos)
K = 4  
result = kmeans(matriz_solar, K)  # Se aplica K-means y se obtiene la matriz con solo 4 casos.

# Obtener los centroides (días típicos)
centroides = result.centers  

# Verificar las dimensiones de la matriz reducida
println("Dimensiones de centroides: ", size(centroides))

# Crear un DataFrame con los centroides (opcional, para visualización)
dias_tipicos = DataFrame(centroides, Symbol.("DT " .* string.(1:K)))

# Mostrar la tabla de centroides
pretty_table(dias_tipicos)

# Se guarda el dataframe de dias_tipicos en un CSV
CSV.write("dias_tipicos.csv", dias_tipicos)

# Se gráfica las curvas de generación por día para los 365 días.
matriz_casos = Matrix(dias_tipicos)
p2 = plot(size=(800, 600), legend=false)
for i in 1:K
    plot!(p2, matriz_casos[:, i], label = false, alpha=0.5)
end
title!("Potencia Generada por día para solo 4 casos típicos")
xlabel!("Intervalos de 5 min")
ylabel!("Potencia Generada")
p3 = plot!(p1,p2, layout = (2,1))
display(p3)

# Se ajustan las demandas a valores promedio cada 5 minutos.

# Promediando los perfiles de demanda que estan por minuto (agrupando cada 5 minutos)
Demandas_df.grupo = div.(0:(nrow(Demandas_df)-1), 5) .+1
Demandas_prom = combine(groupby(Demandas_df, :grupo), names(Demandas_df, Not([:grupo])) .=> mean .=> names(Demandas_df, Not([:grupo])))
select!(Demandas_prom, Not(1))

# Se renombran las columnas de demandas _prom para que coincidan con los nodos PQ.
encabezados = String[]

for i in 1:nrow(nodes)
    if nodes.PLOAD[i] != 0
        push!(encabezados, "N$i")
    end
end

matriz_dem_prom = Matrix(Demandas_prom)
Demandas_prom = DataFrame(matriz_dem_prom, encabezados)

# Se guarda el dataframe de Demanda_prom en un CSV
CSV.write("Demandas_5min.csv", Demandas_prom)

# Grafica de Demanda promedio de cada nodo en intervalos de 5 min.
pbase = 100
p4 = plot(size=(800, 600), legend=false)
for i in 1:(ncol(Demandas_df)-1)
    plot!(p4, matriz_dem_prom[:, i]*pbase, label = false, alpha=0.5)
end
title!("Potencia Demandada por día")
xlabel!("Intervalos de 5 min")
ylabel!("Potencia Demandada")
display(p4)


## FLUJO CUASI-DINAMICO


function Flujo_cuasi_din(lines, nodes, Demandas_prom, dias_tipicos)
    """
    Entradas:   lines: DataFrame con los datos de las lineas del sistema
                nodes : DataFrame con los datos de los nodos del sistema
    Salida :    Ninguna. Nota: La función crea K archivos csv con las tensiones nodales para cada día típico.

    """
    num_nodes = nrow(nodes)
    # Cálculo de Ynn y Yns
    Ynn, Yns = calcular_Ynn_Yns(lines,nodes)
    iYnn = inv(Ynn)
    
    # Se calcula el flujo cuasi-dinámico
    for i in 1:ncol(dias_tipicos)   
        Vns = DataFrame()                 # Se itera en el número de casos
        for j in 1:nrow(dias_tipicos)                # Se itera en el número de datos por intervalos de 5 min
            # Flujo Estatico (intervalo de 5 min)
            # Se inicializa Vn y Vs.
            Vn = ones(num_nodes-1) + 1im*zeros(num_nodes-1)
            Vs = 1 + 0*1im
            Sn = zeros(num_nodes-1)
            for k in 1:ncol(Demandas_prom)           # Se itera en el número de nodos que contienen demandas.
                nodo = parse(Int64,names(Demandas_prom)[k][2:end]) -1
                if  (nodo + 1) == 34                                         # Se escoje el nodo 34 como nodo de generación.
                    Sn[nodo] = dias_tipicos[j,i] - Demandas_prom[j,k]
                end
                Sn[nodo] = - Demandas_prom[j,k]
            end
            # Se inicia el punto fijo.
            for i= 1:4
                # Ecuación de Vn para punto fijo.
                Vn = iYnn*(conj.(Sn./Vn) .- Yns*Vs) 
            end
            # Se guardan los valores en un Dataframe
            nombre_col = "int $(SolarData_df.Hora[j])"
            Vns[!, nombre_col[1:9]] = Vn
        end
        nombre_archivo = "Vn_caso_$i".*".csv"
        CSV.write(nombre_archivo, Vns)
    end
end

# Se ejecuta el flujo cuasi dinamico
Flujo_cuasi_din(lines, nodes, Demandas_prom, dias_tipicos)