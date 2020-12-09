using JLD2
using FileIO

using Glob
using LinearAlgebra


include("simulations-bcel-DM-diss.jl")


function give_alpha_set(n,αmax)
    return rand(n).*αmax
end


function takememory_alpha_set(β,v,α)

    prob=probmemory_alpha_set(β,v,α)
    wahrscheinlikeit = rand()
    result=false
    if wahrscheinlikeit < prob
        result=true
    end
    result
end


function probmemory_alpha_set(β,v,α)
    prob = 1- exp(-β *mean(energymemory.(v,α)))
    prob
end

function onestep_dritf_real_alpha_set(β,δ,oldv,αs,n)
    pos=0
    if (takememory_alpha_set(β,oldv,αs))
        probs=energymemory.(oldv,αs)
        d = Categorical(probs./sum(probs))
        pos=rand(d)
        α=αs[pos]
        free=freeenergy(β,oldv,α)
        diss=dissipation(β,oldv,α)
        unet=free -diss
        oldv=newv_drift(δ,n)
        n+=1

    else
        unet = 0
        free=0
        diss=0
        oldv=newv_drift(δ,1)

        n=2
    end

    return (unet,oldv,free,diss,n,pos)
end


function livesycle_drift_real_alpha_set(β,δ,αs,meanencounters,reps)
    meanunet=0 #initializing for all the rounds for fist cost
    meanfree=0
    meandiss=0
    naiveuses = 0
    for round in 1:reps
        encounters = meanencounters *1
        localnetu= 0
        netu = 0
        localfree= 0
        free = 0
        localdiss= 0
        diss = 0


        v=newv_drift(δ,1) #start with memory produced at t= 0
        naiveuses += 1/reps/encounters
        n=2
        for encounter in 2:encounters  #start at 1 because initialization is allways at 0 in the first round without any prior memory
             naiveuse =0    # Used naive for the first virus
            (netu,v,free,diss,n,dummi)=onestep_dritf_real_alpha_set(β,δ,v,αs,n)
            localnetu+=netu
            localfree+=free
            localdiss+=diss
        end

        meanunet+= localnetu/reps/encounters
        meanfree+= localfree/reps/encounters
        meandiss+= localdiss/reps/encounters
    end

    return (meanunet*meanencounters,meanfree*meanencounters,meandiss*meanencounters)
end


function livesycle_drift_real_alpha_set_deltauniform(β,δmax,αs,meanencounters,reps)
    meanunet=0 #initializing for all the rounds for fist cost
    meanfree=0
    meandiss=0
    for round in 1:reps
        δ=rand()*δmax
    localnetu,localfree,localdiss=livesycle_drift_real_alpha_set(β,δ,αs,meanencounters,1)



        meanunet+= localnetu/reps
        meanfree+= localfree/reps
        meandiss+= localdiss/reps
    end

    return (meanunet,meanfree,meandiss)
end

function alpharange(α,amin,amax) # puts alpha in to the range between alpha_min and alpha_max
    return max(amin,min(α,amax))
end
function optimization_set_alphas(β,δmax,αs,meanencounters,reps,steps,gradients,λ1,λ2)  ## does "steps" optimization steps starting from αs; for  αs==zeros(l) gives rand alpha of langth(l)
    amax=4
    amin=0.001
    l=length(αs)
    if (αs[1]==0)
        αs=rand(l)*amax
    end
    αs=alpharange.(αs,amin,amax)
    astart=αs.*1
    Estart,dunmi1,dumi2=livesycle_drift_real_alpha_set_deltauniform(β,δmax,αs,meanencounters,reps)
    for ii in 1:steps
        E0,dunmi1,dumi2=livesycle_drift_real_alpha_set_deltauniform(β,δmax,αs,meanencounters,reps)
        gradient=zeros(l)
        for jj in 1:gradients
            alocal=alpharange.(αs.+(λ1*randn(l)),amin,amax)
            E1,dunmi1,dumi2=livesycle_drift_real_alpha_set_deltauniform(β,δmax,alocal,meanencounters,reps)
            if (norm((alocal.-αs),2) >0)
                gradient=gradient.+(E1-E0)*(alocal.-αs)/norm((alocal.-αs),2)/gradients
            end
        end
        αs=alpharange.(αs.+λ2*gradient,amin,amax)
    end
    Efinal,dunmi1,dumi2=livesycle_drift_real_alpha_set_deltauniform(β,δmax,αs,meanencounters,reps)
    return αs,astart,Estart,Efinal
end


function simulation_alphaopt_rep(β,δmax,alphas,meanencounters,reps,steps,gradients,λ,λ2,rep)




  ##αs==zeros(l) gives rand alpha of langth(l)
parameters=(β,δmax,alphas,meanencounters,reps,steps,gradients,λ,λ2,rep)

    filename = "Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)/reps_$(reps)/beta_$(β)-Delmax_$(δmax)-asteps_$(alphas)-grad_$(gradients)_steps_$(steps)_l1_$(λ)_l2_$(λ2)_rep_$(rep).jld2"
    println("I want to do  $(filename)")
    if isfile(filename)
        println("Die daten $(filename) exestieren schon ich hoehre auf")
    else

        if !isdir("Data/DM_bcell_drift_alphaopt")
            mkdir("Data/DM_bcell_drift_alphaopt")
        end

        if !isdir("Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)")
            mkdir("Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)")
        end
        if !isdir("Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)/reps_$(reps)")
            mkdir("Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)/reps_$(reps)")
        end

        afin,astat,Estart,efin=optimization_set_alphas(β,δmax,zeros(alphas),meanencounters,reps,steps,gradients,λ,λ2)


        save(filename,"Astart",  astat,"Afin",afin,"Estart",Estart,"Efin",efin,"parameters",parameters)
        println("I'm saving the file $(filename)")
    end
end

function simulation_alphaopt_rep_nextsteps(β,δmax,alphas,meanencounters,reps,steps,gradients,λ,λ2,rep,startstep)




  ##αs==zeros(l) gives rand alpha of langth(l)
parameters=(β,δmax,alphas,meanencounters,reps,steps,gradients,λ,λ2,rep)

    filename = "Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)/reps_$(reps)/beta_$(β)-Delmax_$(δmax)-asteps_$(alphas)-grad_$(gradients)_steps_$(startstep)_l1_$(λ)_l2_$(λ2)_rep_$(rep).jld2"
    if isfile(filename)

        astart_original= load(filename,"Astart")
        afin_original= load(filename,"Afin")

        println("Loaded data and will start simulation")

        filename = "Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)/reps_$(reps)/beta_$(β)-Delmax_$(δmax)-asteps_$(alphas)-grad_$(gradients)_steps_$(steps)_l1_$(λ)_l2_$(λ2)_rep_$(rep).jld2"
        println("I want to do  $(filename)")
        if isfile(filename)
            println("Die daten $(filename) exestieren schon ich hoehre auf")
        else

            if !isdir("Data/DM_bcell_drift_alphaopt")
                mkdir("Data/DM_bcell_drift_alphaopt")
            end

            if !isdir("Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)")
                mkdir("Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)")
            end
            if !isdir("Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)/reps_$(reps)")
                mkdir("Data/DM_bcell_drift_alphaopt/meanencounters_$(meanencounters)/reps_$(reps)")
            end
            rounds=round(Int64,steps-startstep)
            #println("ich starte mit $((β,δmax,afin_original,meanencounters,reps,rounds,gradients,λ,λ2))")
            afin,astat,Estart,efin=optimization_set_alphas(β,δmax,afin_original,meanencounters,reps,rounds,gradients,λ,λ2)
            #println("ich bin fertig mit mit $((β,δmax,afin_original,meanencounters,reps,rounds,gradients,λ,λ2))")
        #filename = "/Users/oskar/Dropbox/PHD/Simulations/b-cell/Simulation/src/julia/data/test.jld"
            save(filename,"Astart",  astart_original,"Afin",afin,"Estart",Estart,"Efin",efin,"parameters",parameters)
            println("I'm saving the file $(filename)")
        end


    else
        println("Data $(filename) allready exists")
    end

end

# Single alpha
