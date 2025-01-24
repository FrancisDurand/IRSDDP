using Random
using Distributions 
using Plots
using JuMP
using CPLEX
using LinearAlgebra


function buildBh_1d!(Bh)
    Emax = 300 #stock max IL FAUDRA CHANGER
    MOI.set(Bh, MOI.Silent(), true)

    @variable(Bh, 0 <= ψ)#Variable epigraphique de Bh
    @variable(Bh, x[1:1])

    #Définir X notre espace d'état
    @constraint(Bh, 0 <= x[1])
    @constraint(Bh, x[1] <= Emax)

    @objective(Bh, Min, ψ)

    Bh.ext[:cuts] = []
    Bh.ext[:ψ] = ψ # La variable epigraphique
    Bh.ext[:x] = x
    Bh.ext[:prec_argmin] = 0 #stoque le dernier argmin calculé de Bh
end

function build_1d!(m)
    nb_palier_th = 4 #nombre de paliers pour la production d'énergie thermique

    #en unité d'énergie
    d = 80 #Demande sur un pas de temps
    gh_sup = 100 #Génération hydrolique max sur un pas de temps
    gt_sup = 40 #Génération thermique max sur un pas de temps
    Emax = 300 #stock max
    capa_th = [10 15 10 5]

    #En coût
    K = 50 #coût d'une unité de blackout
    c_th = [1 2 5 10]
    c_spill = 1
    B_inf = 0 # - 1e8 #Borne inf, en principe il ne faut pas qu'elle soit atteinte
    # x_et = 150

    beta = 0.9

    x_0 = 100
    val_arrivee = 50
    coupes = []
    saa = [50 90]
    nb_mc = length(saa)
    MOI.set(m, MOI.Silent(), true)

    @variable(m, 0 <= z[1:1] <= Emax)
    # @variable(m, z) # position à l'avenir
    # @constraint(m, X, 0 <= z <= Emax) # Domaine de définition de l'avenir

    @variable(m, 0 <= gh[1:1] <= gh_sup)
    # @variable(m, 0 <= gt[1:1][i=1:nb_palier_th] <= capa_th[i])
    @variable(m, 0 <= gt[1, i=1:nb_palier_th] <= capa_th[i])
    @variable(m, 0 <= gt_tot[1:1])
    @variable(m, 0 <= spill[1:1]) #quantité d'eau relachée sans être turbinée
    @variable(m, 0 <= def[1:1]) #demande non satisfaite
    @variable(m, arrivee[1:1])
    @variable(m, x[1:1])

    #Variable pour encoder la suite : V(z)
    @variable(m, Vz >= B_inf) #Je borne la V(z) par 0 pour que V soit définie partout mais c'est pas clair pour moi que c'est compatible avec RVI

    #contrainte pour les coupes
    # @constraint(m, cst_coupe, x - x_0 == 0)
    fix(x[1],x_0)
    fix(arrivee[1], val_arrivee)


    #Contraintes du problème
    # @constraint(m, c_tot, sum(gt[1][i] for i in 1:nb_palier_th) == gt_tot[1])
    @constraint(m, c_tot, sum(gt[1, i] for i in 1:nb_palier_th) == gt_tot[1])
    for i in 1:1
        @constraint(m, cz, z[i] == x[i] - gh[i] + arrivee[i] - spill[i])
        @constraint(m, cd, d[i] == gh[i] + gt_tot[i] + def[i])
    end
    
    #Contraintes des coupes 
    for (x_c,y_c,v_c) in coupes
        @constraint(m, Vz >= y_c + dot(v_c,(z-x_c))) #Ici il faudra remplacer par les mêmes opérations mais vectorielles
    end

    @expression(m, imm_cost, sum(gt[1, i]*c_th[i] for i in 1:nb_palier_th) + K*def[1] + spill[1]*c_spill) #changement

    @objective(m, Min, imm_cost + beta*Vz)

    m.ext[:cuts] = []
    m.ext[:Θ] = Vz # La variable epigraphique
    m.ext[:Β] = beta
    m.ext[:x] = x
    m.ext[:z] = z
    m.ext[:B_inf] = B_inf
    m.ext[:arrivee] = arrivee
    m.ext[:saa] = saa
    m.ext[:x_0] = x_0
    m.ext[:imm_cost] = imm_cost
    m.ext[:dimension_X] = 1
    m.ext[:period] = 1
    

    # m.ext[:X] = X #L'ensemble d'états

    Bh = JuMP.Model(CPLEX.Optimizer)
    buildBh!(Bh)
    m.ext[:Bh] = Bh
    optimize!(m)
end
