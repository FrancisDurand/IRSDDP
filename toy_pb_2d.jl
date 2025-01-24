using Random
using Distributions 
using Plots
using JuMP
using CPLEX
using LinearAlgebra


function buildBh_2d!(Bh, Emax)
#     Emax_1 = 300 #stock max IL FAUDRA CHANGER
#     Emax_2 = 500 #stock max IL FAUDRA CHANGER
    MOI.set(Bh, MOI.Silent(), true)

    @variable(Bh, 0 <= ψ)#Variable epigraphique de Bh
    @variable(Bh, x[1:2])

    #Définir X notre espace d'état
    @constraint(Bh, 0 <= x[1:2])
    @constraint(Bh, x[1] <= Emax[1])
    @constraint(Bh, x[2] <= Emax[2])

    @objective(Bh, Min, ψ)

    Bh.ext[:cuts] = []
    Bh.ext[:ψ] = ψ # La variable epigraphique
    Bh.ext[:x] = x
    Bh.ext[:prec_argmin] = [0, 0] #stoque le dernier argmin calculé de Bh
end

# trans_max is the maximum quantity of energie that can go from dam 1 to dam 2, the cost is trans_c
function build_2d!(m, nb_palier_th,d, gh_sup, Emax, K, c_spill, alpha, x_0, val_arrivee, coupes, saa, capa_th, c_th, trans_max, trans_c)
    MOI.set(m, MOI.Silent(), true)
    nb_mc = length(saa[1])
    B_inf = 0 # - 1e8 #Borne inf, en principe il ne faut pas qu'elle soit atteinte
    # @variable(model, x[1:2] <= ub)

    # @variable(model, x[i=1:2] <= ub[i])

    @variable(m, 0 <= z[j=1:2] <= Emax[j])
    # @variable(m, z) # position à l'avenir
    # @constraint(m, X, 0 <= z <= Emax) # Domaine de définition de l'avenir

    @variable(m, 0 <= gh[j=1:2] <= gh_sup[j])
    @variable(m, 0 <= trans[j=1:2] <= trans_max[j])
    #trans[j] -> quantity going from j to 3-j
    
    # println(capa_th)
    @variable(m, 0 <= gt[j = 1:2, i = 1:nb_palier_th[j]] <= capa_th[j][i])
    # @variables(m, gt[j = 1:2, i = 1:nb_palier_th[j]])
    
    # for j in 1:2
    #     for i in 1:nb_palier_th[j]
    #         set_lower_bound(gt[j][i],0)
    #         println(capa_th[j][i])
    #         set_upper_bound(gt[j][i], capa_th[j][i])
    #     end
    # end

    # @variable(m, 0 <= gt[1, i=1:nb_palier_th[1]] <= capa_th_1[i])
    # @variable(m, 0 <= gt[2, i=1:nb_palier_th[2]] <= capa_th_2[i])
    @variable(m, 0 <= gt_tot[j=1:2])
    @variable(m, 0 <= spill[j=1:2]) #quantité d'eau relachée sans être turbinée
    @variable(m, 0 <= def[1:2]) #demande non satisfaite
    @variable(m, arrivee[1:2])
    @variable(m, x[1:2])
    

    #Variable pour encoder la suite : V(z)
    @variable(m, Vz >= B_inf) #Je borne la V(z) par 0 pour que V soit définie partout mais c'est pas clair pour moi que c'est compatible avec RVI

    #contrainte pour les coupes
    # @constraint(m, cst_coupe, x - x_0 == 0)

    for j in 1:2
        fix(x[j], x_0[j])
        # println(FixRef.(x[j]))
    end
    # fix(x[j = 1:2], x_0[j])
    for j in 1:2
        fix(arrivee[j], val_arrivee[j])
    end
    # fix(arrivee[1:2], val_arrivee)


    #Contraintes du problème
    # @constraint(m, c_tot, sum(gt[1][i] for i in 1:nb_palier_th) == gt_tot[1])
    @constraint(m, c_tot[j in 1:2], sum(gt[j, i] for i in 1:nb_palier_th[j]) == gt_tot[j])

    # for j in 1:2
    #     @constraint(m, c_tot, sum(gt[j][i] for i in 1:(nb_palier_th[j])) == gt_tot[j])
    # end
    # @constraint(m, c_tot, sum(gt[1, i] for i in 1:nb_palier_th) == gt_tot[1])
    for j in 1:2
        @constraint(m, z[j] == x[j] - gh[j] + arrivee[j] - spill[j] - trans[j] + trans[3-j])
        @constraint(m, d[j] == gh[j] + gt_tot[j] + def[j])
    end
    
    #Contraintes des coupes 
    for (x_c,y_c,v_c) in coupes
        @constraint(m, Vz >= y_c + dot(v_c,(z-x_c)))
    end

    @expression(m, imm_cost, sum(sum(gt[j, i]*c_th[j][i] for i in 1:nb_palier_th[j]) + K[j]*def[j] + spill[j]*c_spill[j] + trans[j]*trans_c for j in 1:2))

    @objective(m, Min, imm_cost + alpha*Vz) 

    m.ext[:cuts] = []
    m.ext[:Θ] = Vz # La variable epigraphique
    m.ext[:Β] = alpha
    m.ext[:x] = x
    m.ext[:z] = z
    m.ext[:B_inf] = B_inf
    m.ext[:arrivee] = arrivee
    m.ext[:saa] = saa
    m.ext[:x_0] = x_0
    m.ext[:imm_cost] = imm_cost
    m.ext[:dimension_X] = 2
    # m.ext[:period] = 1
    m.ext[:nb_mc] = nb_mc
    # m.ext[:X] = X #L'ensemble d'états

    Bh = JuMP.Model(CPLEX.Optimizer)
    buildBh_2d!(Bh, Emax)
    m.ext[:Bh] = Bh
    optimize!(m)
end

#This functions builds $period$ models,, list_d and list_saa are the lists of demand and random events that occur.
function build_2d_periodic(period, nb_palier_th, list_d, gh_sup, Emax, K, c_spill, alpha, x_0, val_arrivee, list_saa, capa_th, c_th, trans_max, trans_c)
    models = [JuMP.Model(CPLEX.Optimizer) for _ in 1:period]
    for i in 1:period
        build_2d!(models[i], nb_palier_th, list_d[i], gh_sup, Emax, K, c_spill, alpha, x_0, val_arrivee, [], list_saa[i], capa_th, c_th, trans_max, trans_c)
    end
    return(models)
end