using Poskus

using SparseArrays
using LightGraphs
using Plots

"""
Najprej se izračuna matrika povezanosti A iz grafa g. 
Nato se ta matrika A preoblikuje v sistem enačb, 
kjer se vsaka vrstica prilagodi glede na število sosednjih vozlišč.

Naključno se določijo začetne vrednosti x in y koordinat 
za vsako vozlišče.
Nato se za vozlišča z liho številko izračunajo nove x in y 
koordinate z uporabo SOR metode, medtem ko so za vozlišča 
z sodim indeksom koordinate nastavljene na nič.

Funkcija vrne končne x in y koordinate ter število iteracij, 
potrebnih za dosego konvergence.

"""

function graf(g, omega)
    A2 = adjacency_matrix(g)
    A = copy(A2)

    br = 1
    for i in 1:size(A, 1)
        if isodd(br)
            suma = sum(A[i, :])
            A[i, i] = -suma
        else
            A[i, :] .= 0
            A[i, i] = 1
        end
        br += 1
    end
    
    n = nv(g)
    
    x_zacetni = rand(nv(g)) .- 0.5
    y_zacetni = rand(nv(g)) .- 0.5

    bx = zeros(size(A, 1))
    by = zeros(size(A, 1))
    br = 1
    for i in 1:size(A, 1)
        if isodd(br)
            bx[i] = x_zacetni[i]
            by[i] = y_zacetni[i]
        end
        br += 1
    end
    
    x, iter = sor(RazprsenaMatrika(A), bx, x_zacetni, omega)
    y, iter = sor(RazprsenaMatrika(A), by, y_zacetni, omega)
    
    return x, y, iter
end





g = SimpleGraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 1)
add_edge!(g, 1, 3)
add_edge!(g, 1, 5)
add_edge!(g, 4, 5)


x, y, iter = graf(g, 2)

omega = [0.5, 1, 1.1, 1.5, 1.7, 2]
iter = [62, 16, 11, 37, 78, 500]

# Vizualizacija odvisnosti števila iteracij od omega
plot(omega, iter, xlabel="omega", ylabel="broj iteracij", title="Odvisnost konvergence hitrosti glede na omega.")

# Visualizacija grafa
scatter(x, y, title="Graf v prostoru", markersize=4)
for e in edges(g)
    i, j = src(e), dst(e)
    plot!([x[i], x[j]], [y[i], y[j]], color=:black, linewidth=2)
end
display(plot!(legend=false, xlabel="X", ylabel="Y"))