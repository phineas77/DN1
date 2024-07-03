module Poskus

export RazprsenaMatrika
export sor

import Base.getindex
import Base.setindex!
import Base.firstindex
import Base.lastindex
import Base.*

using LinearAlgebra

"""
`RazprsenaMatrika(V, I)` je poseben podatkovni tip, ki se uporablja za predstavitev razpršenih matrik.
V matriki 'V' so shranjeni neničelni elementi izvirne matrike, 
medtem ko matrika 'I' hrani indekse stolpcev teh neničelnih elementov. 
Vsaka vrstica matrike 'V' vsebuje neničelne elemente iz iste vrstice kot izvirna matrika. 
Če je izvirna matrika dimenzije n x n, potem sta dimenziji 'V' in 'I' enaki n x m.
"""

# Definicija podatkovnega tipa RazprsenaMatrika
mutable struct RazprsenaMatrika
    V  # Matrika n x m, ki vsebuje neničelne elemente matrike A
    I::Matrix{Int64}  # Matrika n x m, ki vsebuje indekse stolpcev neničelnih elementov
end

"""
`RazpršenaMatrika(matrika)` je funkcija, ki sprejme diagonalno dominantno razpršeno matriko
kot vhod in jo pretvori v podatkovni tip RazprsenaMatrika. 
Izhod te funkcije sta matriki 'V' in 'I', ki vsebujeta neničelne elemente in 
njihove ustrezne indekse stolpcev.

Na primer:
    Za vhodno matriko A =[3 0 -2 1 0]
                         [0 3 1 2 0]
                         [0 0 0 0 0]
                         [-1 0 0 2 0]
                         [0 4 0 0 7],
    bosta izhodni matriki:
        V: [3 -2 1]
           [3 1 2]
           [0 0 0]
           [-1 2 0]
           [4 7 0]
        
        I: |1 3 4|
           |2 3 4|
           |0 0 0|
           |1 4 0|
           |2 5 0|

"""

function RazprsenaMatrika(matrika)
    n = size(matrika, 1)
    m = size(matrika, 2)

    # Ustvari matriko I z indeksi stolpcev neničelnih elementov
    I = zeros(Int, n, m)

    # Napolni matriko I z indeksi neničelnih elementov
    for i in 1:n
        k = 1  # Kazalec za shranjevanje indeksov stolpcev neničelnih elementov
        for j in 1:m
            if !iszero(matrika[i, j])
                I[i, k] = j
                k += 1
            end
        end
    end

    # Pridobi maksimalno število neničelnih elementov v vrsticah
    max_non_zero_per_row = maximum([sum(!iszero, matrika[i, :]) for i in 1:n if sum(!iszero, matrika[i, :]) > 0])

    # Ustvari matriko V z ustrezno velikostjo
    V = zeros(n, max_non_zero_per_row)

    # Napolni matriko V z neničelnimi elementi iz matrike A
    for i in 1:n
        for j in 1:max_non_zero_per_row
            idx = I[i, j]
            if idx != 0
                V[i, j] = matrika[i, idx]
            else
                break
            end
        end
    end

    # Ustvari matriko I z ustrezno velikostjo
    I_trimmed = zeros(Int, n, maximum([findlast(!iszero, I[i, :]) for i in 1:n if sum(!iszero, I[i, :]) > 0]))

    # Izreži I, da odstranimo morebitne ničle na koncu vsake vrstice
    for i in 1:n
        last_nonzero_idx = findlast(!iszero, I[i, :])
        if last_nonzero_idx != nothing
            I_trimmed[i, 1:last_nonzero_idx] = I[i, 1:last_nonzero_idx]
        end
    end

    return RazprsenaMatrika(V, I_trimmed)
end

"""
Funkcija `getindex(A::RazprsenaMatrika, i::Int, j::Int)` sprejme tri vhodne parametre: matriko tipa RazprsenaMatrika, 
indeks vrstice in indeks stolpca. Nato vrne element matrike, ki se nahaja na določeni poziciji,
določeni z indeksom vrstice in stolpca.
"""

function getindex(A::RazprsenaMatrika, i::Int, j::Int)
    # Preveri, ali sta indeksa i in j v mejah veljavnosti matrike
    if i < 1 || i > size(A.V, 1) || j < 1 || j > size(A.V, 1)
        throw(ArgumentError("Invalid index"))
    end
    
    # Poišči ustrezen element glede na strukturo razpršene matrike
    for k = 1:size(A.I, 2)
        if A.I[i, k] == j 
            return A.V[i,k]
        end
    end
    
    # Če elementa ni mogoče najti, vrnemo privzeto vrednost (običajno 0 za razpršene matrike)
    return 0
end

"""
Funkcija `firstindex(A::RazprsenaMatrika)` poišče in vrne prvi indeks v matriki, ki je tipa RazprsenaMatrika.
"""

function firstindex(A::RazprsenaMatrika)
    return 1, 1
end


"""
Funkcija lastindex(A::RazprsenaMatrika) vrne zadnji indeks matrike, ki je tipa RazprsenaMatrika. 
"""

# lastindex metoda za RazprsenaMatrika
function lastindex(A::RazprsenaMatrika)
    return size(A.V, 1), size(A.V, 1)
end

"""
Funkcija `*(A::RazprsenaMatrika, x::Vector)` sprejme matriko tipa RazprsenaMatrika in vektor ter izvede njuno matrično množenje. 
"""

# "*" operator za multiplikacijo RazprseneMatrike z vektorjem
function *(A::RazprsenaMatrika, x::Vector)
    if size(A.V, 1) != length(x) 
        throw(error("Velikost matrike in vektorja nista združljiva!"))
    end

    b = zeros(length(x))
    for i=1:size(A.V, 1)
        for j=1:size(A.V, 2)
            if A.I[i, j] != 0
                b[i] = b[i] + A.V[i, j]*x[A.I[i, j]]
            end
        end
    end
    return b
end


"""
Funkcija `sor(A::RazprsenaMatrika, b::Vector, x0::Vector, omega, tol=1e-10)` izvede reševanje
sistema linearnih enačb Ax=b z uporabo metode Successive-over-relaxation (SOR). Vhodni parametri funkcije so:

A: matrika tipa RazprsenaMatrika,
b: vektor, ki predstavlja desno stran sistema enačb,
x0: začetni približek rešitve,
omega: relaksacijski parameter,
tol: toleranca, ki določa pogoj za končanje iteracij.

Funkcija iterativno izvaja SOR iteracije, dokler ni dosežena želena toleranca. 
Na koncu vrne vektor x, ki predstavlja približno rešitev sistema enačb Ax=b.
"""

function sor(A::RazprsenaMatrika, b::Vector, x0::Vector, omega, tol=1e-10)
    max_it = 500
    for it=1:max_it
        for i=1:size(A.V, 1)
            suma = 0
            for k=1:size(A.I, 2)
                if i != A.I[i, k] && A.I[i, k] != 0
                    suma = suma + A.V[i, k]*x0[A.I[i, k]]
                end
            end
            x0[i] = (1 - omega) * x0[i] + omega / getindex(A, i, i) * (b[i] - suma)
        end
        if norm(A*x0 - b) < tol
            return x0, it
        end
    end
    return x0, max_it
    
end


"""
Funkcija `setindex!(A::RazprsenaMatrika, vrednost, i::Int, j::Int)` sprejme matriko tipa RazprsenaMatrika, 
želeno vrednost ter indeks vrstice in indeks stolpca, na katerem želimo shraniti vrednost. 
Funkcija nato vstavi podano vrednost na določeno pozicijo v matriki.
"""

function setindex!(A::RazprsenaMatrika, vrednost, i::Int, j::Int)
    if i > size(A.V, 1) || j > size(A.V, 1) || i < 1 || j < 1
        throw(BoundsError(A, (i, j)))
    end
    
    for k = 1:size(A.I, 2)
        if A.I[i, k] == j
            A.V[i, k] = vrednost
            return A
        end
    end
    
    if j < A.I[i, 1] || (A.I[i, size(A.I, 2)] != 0 && j > A.I[i, end])
        tempV = zeros(size(A.V, 1), size(A.V, 2) + 1)
        tempI = zeros(size(A.V, 1), size(A.V, 2) + 1)
        for k = 1:size(tempV, 1)
            for n = 1:size(tempV, 2)
                if k == i && n == 1
                    tempV[k, n] = vrednost
                    tempI[k, n] = j
                elseif k == i && n == size(tempV, 2)
                    tempV[k, n] = A.V[k, end]
                    tempI[k, n] = A.I[k, end]
                elseif k == i
                    tempV[k, n] = A.V[k, n-1]
                    tempI[k, n] = A.I[k, n-1]
                elseif n == size(tempV, 2)
                    tempV[k, n] = 0
                    tempI[k, n] = 0
                else
                    tempV[k, n] = A.V[k, n]
                    tempI[k, n] = A.I[k, n]
                end
            end
        end
        A.V = tempV
        A.I = tempI
        return A
    end
    
    if A.I[i, 1] == 0
        A.I[i, 1] = j
        A.V[i, 1] = vrednost
    else
        for k = 1:size(A.I, 2)
            if A.I[i, k] > j
                A.I[i, k+1:end] = A.I[i, k:end-1]
                A.V[i, k+1:end] = A.V[i, k:end-1]
                A.I[i, k] = j
                A.V[i, k] = vrednost
                break
            end
        end
    end
    return A
end
end # module Poskus
