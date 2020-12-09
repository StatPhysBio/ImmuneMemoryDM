using Distributions
using SpecialFunctions
using JLD2
using FileIO
using Distributed

function minalpha(β)
    α =   (exp(1)* sqrt( π))/β
    return α

end

function takememory(β,v,α)
    prob =Probability_of_memory_usage(β,v,α)
    wahrscheinlikeit = rand()
    result=false
    if wahrscheinlikeit < prob
        result=true
    end
    result
end




function Probability_of_memory_usage(β,v,α)
    prob = 1- exp(-β *energymemory(v,α))
    prob
end

function energymemory(v,α)
    energy =α*exp(-(α*v)^2)*sqrt(1/ π)
    energy
end



function freeenergy(β,v,α)
    energy = log( exp(β *energymemory(v,α)) -1)/β
    energy
end

function dissipation(β,v,α)
    diss = α*α/β * (v)^2
    diss
end





function newv_drift(δ,n)
    vnew=sqrt(n)δ+randn()*δ*0.05*sqrt(n)
    vnew
end



function onestep(β,δ,oldv,α,n)

    if (takememory(β,oldv,α))
        free=freeenergy(β,oldv,α)
        diss=dissipation(β,oldv,α)
        unet=free -diss
        nextv=newv_drift(δ,n)
        n+=1
        naiveuse=0
    else
        unet = 0
        free=0
        diss=0
        nextv=newv_drift(δ,1)
        naiveuse=1
        n=2
    end

    return (unet,nextv,naiveuse,free,diss,n)
end



function livesycle(β,δ,α,meanencounters,reps)
    meanunet=0 #initializing for all the rounds for fist cost
    meanfree=0
    meandiss=0
    naiveuses = 0
    for round in 1:reps
        localnetu= 0
        netu = 0
        localfree= 0
        free = 0
        localdiss= 0
        diss = 0


        v=newv_drift(δ,1) #start the cirus at the memory wich is allways at 0
        naiveuses += 1/reps
        n=2
        for encounter in 2:meanencounters  #start at 1 because initialization is allways at 0 in the first round without any prior memory
             naiveuse =0    # Used naive for the first virus
            (netu,v,naiveuse,free,diss,n)=onestep(β,δ,v,α,n)
            localnetu+=netu
            localfree+=free
            localdiss+=diss
            naiveuses += naiveuse/reps
        end

        meanunet+= localnetu/reps
        meanfree+= localfree/reps
        meandiss+= localdiss/reps
    end

    return (meanunet,naiveuses,meanfree,meandiss)
end








function findoptalphaoflive_drift_real(β,δ,meanencounters,reps,αsteps,αmax)

    #if (minalpha(β) >αmax)
    #    println("amax is to small for the chosen beta")
    #    return 0,0,0,0,0,0
    #end
    start= 0.01
    α = range(start,stop=αmax,length=αsteps)

    data=livesycle.(β,δ,α,meanencounters,reps)
    usesmaive=zeros(αsteps)
    memuses=zeros(αsteps)
    utility=zeros(αsteps)
    free=zeros(αsteps)
    diss=zeros(αsteps)
    for ii in 1:αsteps
        utility[ii] = data[ii][1]
        usesmaive[ii] = data[ii][2]
        free[ii] = data[ii][3]
        diss[ii] = data[ii][4]
    end
    memuses= meanencounters./usesmaive

    amax= argmax(utility)
    optalpha = α[amax]
    Nusesatmax=memuses[amax]
    argmaxn = argmax(memuses)
    maxNpossible=memuses[argmaxn]
    utiliymax = utility[amax]
    freemax=free[amax]
    dissmax=diss[amax]


    return Nusesatmax,utiliymax,optalpha,maxNpossible,freemax,dissmax
end


function find_optimal_alpha_for_range_of_delta(β,meanencounters,reps,αsteps,αmax,deltasteps)

    δ = range(0.0001,stop=1.5/αmax,length= deltasteps)
    #δ = range(0.325,stop=.4,length= deltasteps)./sqrt(αmax)
    memstepsatmax =zeros(deltasteps)
    utilityopt =zeros(deltasteps)
    alphaopt =zeros(deltasteps)
    posnmax =zeros(deltasteps)
    freemax =zeros(deltasteps)
    dissmax =zeros(deltasteps)
    for i =1:deltasteps
        memstepsatmax[i],utilityopt[i],alphaopt[i],posnmax[i],freemax[i],dissmax[i]=findoptalphaoflive_drift_real(β,δ[i],meanencounters,reps,αsteps,αmax)
    end
    return memstepsatmax,utilityopt,alphaopt,posnmax,δ,freemax,dissmax

end


function find_optimal_alpha_for_range_delta_and_beta(meanencounters,reps,αsteps,αmax,deltasteps,betasteps)

    β = range(0.1,length=betasteps,stop=10.0)


    memstepsatmax= zeros(betasteps,deltasteps)
    utilityopt= zeros(betasteps,deltasteps)
    alphaopt= zeros(betasteps,deltasteps)
    posnmax= zeros(betasteps,deltasteps)
    freemax =zeros(betasteps,deltasteps)
    dissmax =zeros(betasteps,deltasteps)
    δ=zeros(deltasteps)
    Threads.@threads for ii in 1:betasteps      # for parallel execution
    # for ii in 1:betasteps                     # non parallel
    memstepsatmax[ii,:],utilityopt[ii,:],alphaopt[ii,:],posnmax[ii,:],δ,freemax[ii,:],dissmax[ii,:] = find_optimal_alpha_for_range_of_delta(β[ii],meanencounters,reps,αsteps,αmax,deltasteps)
    end
    introduction="This gives the data meory_used_at_the_max  utility_atmax  otimal_alpha posible_maximum_memoryuses delta and beta and free_at_the_max and diss_at_the_max the data is Gamma first and delta second index"
 return memstepsatmax,utilityopt,alphaopt,posnmax,δ,β,freemax,dissmax,introduction

end
