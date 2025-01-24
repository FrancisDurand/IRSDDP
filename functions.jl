include("msppy_hydro_thermal.jl")
#En mémoire juste
# function Th(m, x)

#     if is_fixed(m.ext[:x])
#         unfix(m.ext[:x])
#     end
#     fix(m.ext[:x], x)

#     somme_obj = 0
#     somme_dual_val = 0

#     for arrivee in m.ext[:saa]
#         if is_fixed(m.ext[:arrivee])
#             unfix(m.ext[:arrivee])
#         end
#         # print(m)
#         fix(m.ext[:arrivee], arrivee)
#         optimize!(m)
#         # Check the solver status
#         status = termination_status(m)
#         # Check if the solution is optimal
#         if status == MOI.OPTIMAL
#             # Objective value
#             obj_value = objective_value(m)
#             next_x = value.(m.ext[:z])
#             dual_val = dual(FixRef(m.ext[:x]))
#             # println(dual_val)

#             somme_obj += obj_value
#             somme_dual_val += dual_val
#         else
#             display("The model was not solved to optimality.")
#         end
#     end
#     moyenne_obj = somme_obj/length(m.ext[:saa])
#     moyenne_dual_val = somme_dual_val/length(m.ext[:saa])
#     return(moyenne_obj, moyenne_dual_val)
# end

# Fonction qui renvoie (t - 1) modulo period, avec un résultat compris entre 1 et period
function prec_per(t, period)
    return mod(t - 2, period) + 1
end

function next_per(t, period)
    return mod(t, period) + 1
end


function Th(m, x)
    # println(x)
    
    # println(fixed_x)
    for i in 1:m.ext[:dimension_X]
        # if is_fixed(m.ext[:x][i])
        #     unfix(m.ext[:x][i])
        # end
        # println(x[i])
        # println(m.ext[:x][i])
        fix(m.ext[:x][i], x[i])
    end
    fixed_x = FixRef.(m.ext[:x])

    somme_obj = 0
    somme_dual_val = zeros(m.ext[:dimension_X])

    for scenario in 1:m.ext[:nb_mc]
        for i in 1:m.ext[:dimension_X]
            if is_fixed(m.ext[:arrivee][i])
                unfix(m.ext[:arrivee][i])
            end
            # println(m.ext[:saa][i][scenario])
            fix(m.ext[:arrivee][i], m.ext[:saa][i][scenario])
        end

        optimize!(m)
        # print(solution_summary(m; verbose = true))

        # Check the solver status
        status = termination_status(m)
        # Check if the solution is optimal
        if status == MOI.OPTIMAL
            # Objective value
            obj_value = objective_value(m)
            # next_x[i in 1:m.ext[:dimension_X]] = value.(m.ext[:z][i])
            next_x = value.(m.ext[:z][1:m.ext[:dimension_X]])
            
            # dual_val = dual(FixRef(m.ext[:x][1:m.ext[:dimension_X]]))
            dual_val = dual.(fixed_x)
            # println("obj val : ", obj_value)
            # println("dual val : ", dual_val)

            somme_obj += obj_value
            somme_dual_val += dual_val
        else
            display("The model was not solved to optimality.") # TODO un jour enregistrer les modèles qui plantent
        end
    end
    moyenne_obj = somme_obj/length(m.ext[:saa])
    moyenne_dual_val = somme_dual_val/length(m.ext[:saa])
    # println(moyenne_obj, moyenne_dual_val)
    return(moyenne_obj, moyenne_dual_val)
end

function cut!(m, x_cut) #Adds a cut to its own model, used for the stationnary case
    Thx_cut, dual_val = Th(m, x_cut)
    cut = (x_cut, Thx_cut, dual_val)
    push!(m.ext[:cuts], cut)
    # println(value.(m.ext[:z]))
    # println(dual_val)
    # println(value.(m.ext[:z]))
    # println(x_cut)
    # println(dot(dual_val, value.(m.ext[:z]) - x_cut))
    @constraint(m, m.ext[:Θ] >= Thx_cut + dot(dual_val, m.ext[:z] - x_cut))
end

function cut_period!(models, x_cut, t, period) #Adds the cut to the previous model, used for the periodic case
    Thx_cut, dual_val = Th(models[t], x_cut)
    cut = (x_cut, Thx_cut, dual_val)
    push!(models[prec_per(t, period)].ext[:cuts], cut)
    # println(value.(m.ext[:z]))
    # println(dual_val)
    # println(value.(m.ext[:z]))
    # println(x_cut)
    # println(dot(dual_val, value.(m.ext[:z]) - x_cut))
    @constraint(models[prec_per(t, period)], models[prec_per(t, period)].ext[:Θ] >= Thx_cut + dot(dual_val, models[prec_per(t, period)].ext[:z] - x_cut))
end

function nb_cuts_added(m)
    return(length(m.ext[:cuts]))
end

function nb_cuts_total(models)
    s = 0
    for i in eachindex(models)
        s+= nb_cuts_added(models[i])
    end
    return(s)
end

function valeur(m, z)
    pot_y = m.ext[:B_inf] #everything will be positive so 0 is the lowest possible value for me !!! C'est dangereux
    for (x_c, y_c, g_c) in m.ext[:cuts]
        # println("z : ",z)
        # println("x_c : ", x_c)
        actual_y = y_c + dot(g_c, z - x_c)
        if pot_y <= actual_y
            pot_y = actual_y
        end
    end
    return(pot_y)
end

function valeur_period(models, z, t)
    return(valeur(models[t],z))
end

function cut_decal_lambda!(m, x_cut, λ)
    obj_val, dual_val = Th(m, x_cut)
    cut = (x_cut, obj_val - λ, dual_val)
    push!(m.ext[:cuts], cut)
    @constraint(m, m.ext[:Θ] >= obj_val -λ + dot(dual_val, m.ext[:z] - x_cut))
end

function cut_decal_lambda_period!(models, x_cut, λ, t, period)
    obj_val, dual_val = Th(models[t], x_cut)
    cut = (x_cut, obj_val - λ, dual_val)
    push!(models[prec_per(t, period)].ext[:cuts], cut)
    @constraint(models[prec_per(t, period)], models[prec_per(t, period)].ext[:Θ] >= obj_val -λ + dot(dual_val, models[prec_per(t, period)].ext[:z] - x_cut))
end

function cut_decal_y!(m, x_cut, y_decal)
    hy = valeur(m, y_decal)
    Thy, _ = Th(m, y_decal)
    cut_decal_lambda!(m, x_cut, Thy - hy)
end

function cut_decal_y_period!(models, x_cut, y_decal, t, period) # probablement bugée si on prend le cut decal lambda
    # println(y_decal)
    # println(models)
    # println(t)
    # println(models[t])
    hy = valeur_period(models, y_decal, prec_per(t, period))
    Thy, _ = Th(models[t], y_decal)
    # println(Thy, " ", hy, " ", Thy - hy)
    cut_decal_lambda_period!(models, x_cut, Thy - hy, t, period)
    # cut_decal_lambda!(models[t], x_cut, Thy - hy)
end

function cut_Bh!(m, x_cut, obj_val, dual_val)
    push!((m.ext[:Bh]).ext[:cuts], (x_cut, obj_val, dual_val))
    @constraint(m.ext[:Bh], (m.ext[:Bh]).ext[:ψ] >= obj_val + dot(dual_val, (m.ext[:Bh]).ext[:x] - x_cut))
end

function cut_Bh_period!(models, x_cut, obj_val, dual_val, t) #The cuts for Bh are added at its time stamp
    push!((models[t].ext[:Bh]).ext[:cuts], (x_cut, obj_val, dual_val))
    @constraint(models[t].ext[:Bh], (models[t].ext[:Bh]).ext[:ψ] >= obj_val + dot(dual_val, (models[t].ext[:Bh]).ext[:x] - x_cut))
end

function nb_cuts_Bh(m)
    return(length((m.ext[:Bh]).ext[:cuts]))
end

function nb_cut_Bh_total(models)
    s = 0
    for i in eachindex(models)
        s+=nb_cuts_Bh(models[i])
    end
    return(s)
end

function valeur_Bh(m, z)
    pot_y = m.ext[:B_inf] #everything will be positive so 0 is the lowest possible value for me !!! C'est dangereux
    for (x_c, y_c, g_c) in (m.ext[:Bh]).ext[:cuts]
        actual_y = y_c + dot(g_c, z-x_c)
        if pot_y <= actual_y
            pot_y = actual_y
        end
    end
    return(pot_y)
end

function valeur_Bh_per(models, z, t)
    return(valeur_Bh(models[t], z))
end

function argmin_Th(m, eps) #eps est la précision
    w = (m.ext[:Bh]).ext[:prec_argmin]
    Bhw, dual_val = Th(m, w)
    # println(w, " ", Bhw, " ", dual_val, " ", valeur_Bh(m, w))
    while Bhw - valeur_Bh(m, w) > eps
        cut_Bh!(m, w, Bhw, dual_val)
        # println(w, " ", Bhw, " ", dual_val, " ", valeur_Bh(m, w))
        optimize!(m.ext[:Bh])
        w = value.((m.ext[:Bh]).ext[:x])
        Bhw, dual_val = Th(m, w)
    end
    (m.ext[:Bh]).ext[:prec_argmin] = w
    return(w)
end

function argmin_Th_period(models, eps, t)
    return(argmin_Th(models[t], eps))
end

function cut_argmin_Bh!(m, x_cut, eps)
    y = argmin_Th(m, eps)
    cut_decal_y!(m, x_cut, y)
end

function cut_argmin_Bh_period!(models, x_cut, eps, t)
    y = argmin_Th_period(models, eps, t)
    # println(y)
    cut_decal_y_period!(models, x_cut, y, t, length(models))
end


function one_step(m, x, arrivee)
    for i in 1:m.ext[:dimension_X]
        fix(m.ext[:x][i], x[i])
    end

    for i in 1:m.ext[:dimension_X]
        if is_fixed(m.ext[:arrivee][i])
            unfix(m.ext[:arrivee][i])
        end
        # println(m.ext[:saa][i][scenario])
        fix(m.ext[:arrivee][i], arrivee[i])
    end
    optimize!(m)
    # Check the solver status
    status = termination_status(m)
    # Check if the solution is optimal
    if status == MOI.OPTIMAL
        # obj_value = objective_value(m)
        next_x = value.(m.ext[:z])
        current_cost = value.(m.ext[:imm_cost])
        return(next_x, current_cost)
    else
        display("The model was not solved to optimality.")
    end
end

function one_step_period(models, x, arrivee, t)
    return(one_step(models[t], x, arrivee)) #next_x is a point at period t+1
end

function policy_evalution(m, horizon)
    coeff = 1.0
    total_cost = 0.0
    current_x = m.ext[:x_0]
    for i in 1:horizon
        arrivee = Random.rand(m.ext[:saa]) #Pas sur que ce soit ça
        current_x, current_cost = one_step(m, current_x, arrivee)
        # println("current_x : ", current_x, " arrivee : ", arrivee, " current_cost : ", current_cost)
        total_cost += current_cost*coeff
        coeff *= m.ext[:Β]
    end
    return(total_cost)
end

function policy_evaluation_period(models, horizon) # the Horizon is the number of steps i.e. approxiametly nb_cycle*period
    coeff = 1.0
    total_cost = 0.0
    current_x = models[1].ext[:x_0]
    current_time = 1
    for i in 1:horizon
        arrivee = zeros(models[current_time].ext[:dimension_X])
        for j in 1:models[current_time].ext[:dimension_X]
            # arrivee[j] = Random.rand(Random.rand(models[current_time].ext[:saa])) Wtf je prenais la région du brésil au hasard
            arrivee[j] = Random.rand((models[current_time].ext[:saa])[j]) #verifier que c'est bon
        end
        # arrivee = Random.rand(models[current_time].ext[:saa])
        current_x, current_cost = one_step_period(models, current_x, arrivee, current_time)
        # println("time : ", current_time, " current_cost : ", current_cost, " length(arrivee) : ", length(arrivee), " current_x : ", current_x )
        total_cost += current_cost*coeff
        coeff*= (models[current_time]).ext[:Β]
        current_time = next_per(current_time, length(models))
    end
    return(total_cost)
end

function multiple_evaluation(m, nb_eval, horizon)
    total = 0.0
    sum_square = 0.0
    for i in 1:nb_eval
        act_val = policy_evalution(m, horizon)
        total += act_val
        sum_square += act_val^2
    end
    return(total/nb_eval, sqrt(sum_square/nb_eval - (total/nb_eval)^2))
end

function multiple_evaluation_period(models, nb_eval, horizon) #Renvoie la moyenne et la variance
    total = 0.0
    sum_square = 0.0
    for i in 1:nb_eval
        act_val = policy_evaluation_period(models, horizon)
        total += act_val
        sum_square += act_val^2
    end
    return(total/nb_eval, sum_square/nb_eval - (total/nb_eval)^2)
end


function multiple_evaluation_period_return_all(models, nb_eval, horizon) 
    all_costs = []
    for i in 1:nb_eval
        act_val = policy_evaluation_period(models, horizon)
        push!(all_costs, act_val)
    end
    return(all_costs)
end


function is_active(m, cut)

    for i in 1:m.ext[:dimension_X]
        if is_fixed(m.ext[:x][i])
            unfix(m.ext[:x][i])
        end
    end
    for i in 1:m.ext[:dimension_X]
        if is_fixed(m.ext[:z][i])
            unfix(m.ext[:z][i])
        end
    end

    # if is_fixed(m.ext[:z])
    #     unfix(m.ext[:z])
    # end
    # if is_fixed(m.ext[:x])
    #     unfix(m.ext[:x])
    # end
    x_c, y_c, d_c = cut
    temp_objective = @expression(m, m.ext[:Θ] - y_c - dot(d_c, m.ext[:z] - x_c))
    set_objective_function(m, temp_objective)
    optimize!(m)
    w = objective_value(m)
    # println(value(m.ext[:z]))
    # println(value(m.ext[:Θ]))
    # println(value(+y_c + d_c*(m.ext[:z] - x_c)))
    # println(w)
    set_objective_function(m, m.ext[:imm_cost] + m.ext[:Β]*m.ext[:Θ])
    return(w <= 1e-6)
end

function is_active_period(models, cut, t)
    return(is_active(models[t], cut))
end

function nb_active_cuts(m)
    nb_active = 0
    for cut in m.ext[:cuts]
        if is_active(m, cut)
            nb_active += 1
        end
    end
    return(nb_active)
end

function nb_active_cuts_period(models)
    s=0
    for t in eachindex(models)
        s+= nb_active_cuts(models[t])
    end
    return(s)
end

function coupe_function(x,coupe)
    (x_c, y_c, g_c) = coupe
    return(max(y_c + dot(g_c, x - x_c),0))
end

function non_active_cut_removal(m)
    new_cuts = []
    for cut in m.ext[:cuts]
        if is_active(m, cut)
            push!(new_cuts, cut)
        end
    end
    m.ext[:cuts] = new_cuts
end

function is_active_in_Bh(Bh, cut)
    for i in 1:Bh.ext[:dimension_X]
        if is_fixed(Bh.ext[:x][i])
            unfix(Bh.ext[:x][i])
        end
    end
    x_c, y_c, d_c = cut
    temp_objective = @expression(Bh, Bh.ext[:ψ] - y_c - dot(d_c, Bh.ext[:x] - x_c))
    set_objective_function(Bh, temp_objective)
    optimize!(Bh)
    w = objective_value(Bh)
    set_objective_function(Bh, Bh.ext[:ψ])
    return(w <= 1e-6)
end

function non_active_cut_removal_period(models)
    for t in eachindex(models)
        non_active_cut_removal(models[t])
    end
end

function non_active_cut_removal_Bh(m)
    new_cuts = []
    for cut in (m.ext[:Bh]).ext[:cuts]
        if is_active_in_Bh(m.ext[:Bh], cut)
            push!(new_cuts, cut)
        end
    end
    (m.ext[:Bh]).ext[:cuts] = new_cuts
end

function non_active_cut_removal_period_Bh(models)
    for t in eachindex(models)
        non_active_cut_removal_Bh(models[t])
    end
end

function nb_active_cut_Bh(m)
    nb_active_cuts = 0
    for cut in (m.ext[:Bh]).ext[:cuts]
        if is_active_in_Bh(m.ext[:Bh], cut)
            nb_active_cuts += 1
        end
    end
    return(nb_active_cuts)
end

function nb_active_cut_period_Bh(models)
    list_active = []
    for t in eachindex(models)
        push!(list_active, nb_active_cut_Bh(models[t]))
    end
    return(list_active)
end

function nb_cut_per_period(models)
    vect_of_nb = zeros(length(models))
    for i in eachindex(models)
        vect_of_nb[i] = nb_cuts_added(models[i])
    end
    println(vect_of_nb)
end

function test_global_1run(Β)
    modelsRVI = build_ht12month(Β)
    modelsVI = build_ht12month(Β)

    for i in 1:400
        for t in 1:12
            unif_1 = Uniform(0.0, 200717.6)
            unif_2 = Uniform(0.0, 19617.2)
            unif_3 = Uniform(0.0, 51806.1)
            unif_4 = Uniform(0.0, 12744.9)
            cut_here_1 = rand(unif_1)
            cut_here_2 = rand(unif_2)
            cut_here_3 = rand(unif_3)
            cut_here_4 = rand(unif_4)
            x_cut = [cut_here_1, cut_here_2, cut_here_3, cut_here_4]
            # println([cut_here_1, cut_here_2, cut_here_3, cut_here_4])
            cut_argmin_Bh_period!(modelsRVI, x_cut, 1e-6, 13-t)
            cut_period!(modelsVI, x_cut, 13-t, 12)
        end
        if i%200 == 0
            non_active_cut_removal_period(modelsRVI)
            non_active_cut_removal_period(modelsVI)
    
            non_active_cut_removal_period_Bh(modelsRVI)
        end
    end

    nb_cut_RVI = nb_cuts_total(modelsRVI)
    nb_cut_VI = nb_cuts_total(modelsVI)

    nb_eval = 100
    horizon = 1000 #10000 des fois
    score_RVI, _ = multiple_evaluation_period(modelsRVI, nb_eval, horizon)
    score_VI, _ = multiple_evaluation_period(modelsVI, nb_eval, horizon)
    return(nb_cut_RVI, nb_cut_VI, score_RVI, score_VI)
end

function test_global_runs(Β, nb_runs)
    vect_nb_cut_RVI = []
    vect_nb_cut_VI = []
    vect_score_RVI = []
    vect_score_VI = []
    for i in 1:nb_runs
        nb_cut_RVI, nb_cut_VI, score_RVI, score_VI = test_global_1run(Β)
        push!(vect_nb_cut_RVI, nb_cut_RVI)
        push!(vect_nb_cut_VI, nb_cut_VI)
        push!(vect_score_RVI, score_RVI)
        push!(vect_score_VI, score_VI)
    end
    return(vect_nb_cut_RVI,
    vect_nb_cut_VI,
    vect_score_RVI,
    vect_score_VI)
end

function test_global_1run_scaled(Β, scale_factor)
    modelsRVI = build_ht12month_scaled(Β, scale_factor)
    modelsVI = build_ht12month_scaled(Β, scale_factor)

    for i in 1:400
        for t in 1:12
            unif_1 = Uniform(0.0, 200717.6)
            unif_2 = Uniform(0.0, 19617.2)
            unif_3 = Uniform(0.0, 51806.1)
            unif_4 = Uniform(0.0, 12744.9)
            cut_here_1 = rand(unif_1)
            cut_here_2 = rand(unif_2)
            cut_here_3 = rand(unif_3)
            cut_here_4 = rand(unif_4)
            x_cut = [cut_here_1, cut_here_2, cut_here_3, cut_here_4]
            # println([cut_here_1, cut_here_2, cut_here_3, cut_here_4])
            cut_argmin_Bh_period!(modelsRVI, x_cut, 1e-6, 13-t)
            cut_period!(modelsVI, x_cut, 13-t, 12)
        end
        if i%200 == 0
            non_active_cut_removal_period(modelsRVI)
            non_active_cut_removal_period(modelsVI)
            non_active_cut_removal_period_Bh(modelsRVI)
        end
    end

    nb_cut_RVI = nb_cuts_total(modelsRVI)
    nb_cut_VI = nb_cuts_total(modelsVI)

    nb_eval = 100
    horizon = 1000 #10000 des fois
    score_RVI, _ = multiple_evaluation_period(modelsRVI, nb_eval, horizon)
    score_VI, _ = multiple_evaluation_period(modelsVI, nb_eval, horizon)
    return(nb_cut_RVI, nb_cut_VI, score_RVI, score_VI)
end


function multiple_random_cuts_VI(models, nb_cuts_added, time_of_last_added_cut) #Les coupes sont ajoutées à rebours
    time = time_of_last_added_cut
    for _ in 1:nb_cuts_added
        storedEnergy_ub = models[1].ext[:storedEnergy_ub]
        x_cut = [rand(Uniform(0.0, storedEnergy_ub[i])) for i in 1:length(storedEnergy_ub)]
        cut_period!(models, x_cut, time, length(models))
        time = prec_per(time, length(models))
    end
end

function multiple_random_cuts_RVI(models, nb_cuts_added, time_of_last_added_cut) #Les coupes sont ajoutées à rebours
    time = time_of_last_added_cut
    for _ in 1:nb_cuts_added
        storedEnergy_ub = models[1].ext[:storedEnergy_ub]
        x_cut = [rand(Uniform(0.0, storedEnergy_ub[i])) for i in 1:length(storedEnergy_ub)]
        cut_argmin_Bh_period!(models, x_cut, 1e-6, time)
        time = prec_per(time, length(models))
    end
end

function multiple_trajectory_cuts_VI(models, nb_cuts_added, x_first_cut, starting_time) #Les coupes sont ajoutées à rebours
    time = starting_time
    cuts_added_here = [x_first_cut]

    #Une Forward Phase
    for _ in 1:(nb_cuts_added-1)
        # arrivee = Random.rand(models[time].ext[:saa])
        random_index = rand(1:length(models[1].ext[:saa][1]))
        arrivee = [models[time].ext[:saa][i][random_index] for i in 1:(models[1].ext[:dimension_X])]
        current_x, _ = one_step(models[time], cuts_added_here[end], arrivee)
        push!(cuts_added_here, current_x)
        time = next_per(time, length(models))
    end

    #Une Backward Phase
    for x_cut in reverse(cuts_added_here)
        cut_period!(models, x_cut, time, length(models))
        time = prec_per(time, length(models))
    end
end

function multiple_trajectory_cuts_RVI(models, nb_cuts_added, x_first_cut, starting_time) #Les coupes sont ajoutées à rebours
    time = starting_time
    cuts_added_here = [x_first_cut]

    #Une Forward Phase
    for _ in 1:(nb_cuts_added-1)
        # arrivee = Random.rand(models[time].ext[:saa])  #Il y avait une erreur, in faut appeler models[time].ext[:saa], en fonction de time c'est pas pareil
        random_index = rand(1:length(models[1].ext[:saa][1]))
        arrivee = [models[time].ext[:saa][i][random_index] for i in 1:(models[1].ext[:dimension_X])]
        current_x, _ = one_step(models[time], cuts_added_here[end], arrivee)
        push!(cuts_added_here, current_x)
        time = next_per(time, length(models))
    end

    #Une Backward Phase
    for x_cut in reverse(cuts_added_here)
        cut_argmin_Bh_period!(models, x_cut, 1e-6, time)
        time = prec_per(time, length(models))
    end
end

function generate_inflows_one_scenario(models, start_time, horizon)
    inflows_one_scenario = []
    time = start_time
    for _ in 1:horizon
        random_index = rand(1:length(models[1].ext[:saa][1]))
        arrivee = [models[time].ext[:saa][i][random_index] for i in 1:(models[1].ext[:dimension_X])]
        push!(inflows_one_scenario, arrivee)
        time = next_per(time, length(models))
    end
    return(inflows_one_scenario)
end

function generate_inflows_multiple_scenarios(models, start_time, horizon, nb_scenarios)
    return([generate_inflows_one_scenario(models, start_time, horizon) for _ in 1:nb_scenarios])
end

function evaluate_policy_given_scenario(models, starting_time, inflows, max_time)
    coeff = 1.0
    total_cost = 0.0
    current_x = models[starting_time].ext[:x_0]
    current_time = 1
    for i in eachindex(inflows)
        if i <= max_time
            arrivee = inflows[i]
            current_x, current_cost = one_step_period(models, current_x, arrivee, current_time)
            # println("time : ", current_time, " current_cost : ", current_cost, " length(arrivee) : ", length(arrivee), " current_x : ", current_x )
            total_cost += current_cost*coeff
            coeff*= (models[current_time]).ext[:Β]
            current_time = next_per(current_time, length(models))
        end
    end
    return(total_cost)
end

function multiple_evaluation_period_given_scenarios(models, starting_time, list_scenarios, max_time) #Renvoie la moyenne et la variance
    total = 0.0
    sum_square = 0.0
    nb_eval = 0
    for i in eachindex(list_scenarios)
        nb_eval += 1
        act_val = evaluate_policy_given_scenario(models, starting_time, list_scenarios[i], max_time)
        total += act_val
        sum_square += act_val^2
    end
    return(total/nb_eval, sum_square/nb_eval - (total/nb_eval)^2)
end

function test_single_run_TF_active_cuts(Β, k)
    modelsRVI = build_ht12month(Β)
    modelsVI = build_ht12month(Β)

    for i in 1:k
        multiple_trajectory_cuts_RVI(modelsRVI, 12*i, modelsRVI[1].ext[:x_0], 1)
        multiple_trajectory_cuts_VI(modelsVI, 12*i, modelsVI[1].ext[:x_0], 1)
        
        non_active_cut_removal_period(modelsRVI)
        non_active_cut_removal_period(modelsVI)
        non_active_cut_removal_period_Bh(modelsRVI)
    end

    nb_cut_RVI = nb_cuts_total(modelsRVI)
    nb_cut_VI = nb_cuts_total(modelsVI)

    return(nb_cut_RVI, nb_cut_VI)
end

function test_multiple_runs_TF_active_cuts(Β, k, nb_runs)
    vect_nb_cut_RVI = []
    vect_nb_cut_VI = []
    for _ in 1:nb_runs
        nb_cut_RVI, nb_cut_VI = test_single_run_TF_active_cuts(Β, k)
        push!(vect_nb_cut_RVI, nb_cut_RVI)
        push!(vect_nb_cut_VI, nb_cut_VI)
    end
    return(vect_nb_cut_RVI, vect_nb_cut_VI)
end

function multiple_phaseRVI(models, current_k, nb_phase)
    for i in current_k:(current_k+nb_phase-1)
        multiple_trajectory_cuts_RVI(models, 12*i, models[1].ext[:x_0], 1)
        if (rand() < 0.2)
            non_active_cut_removal_period(models)
            non_active_cut_removal_period_Bh(models)
            println("unactive cuts removed")
        end
    end
end

function multiple_phaseVI(models, current_k, nb_phase)
    for i in current_k:(current_k+nb_phase-1)
        multiple_trajectory_cuts_VI(models, 12*i, models[1].ext[:x_0], 1)
        if (rand() < 0.2)
            non_active_cut_removal_period(models)
            println("unactive cuts removed")
        end
    end
end