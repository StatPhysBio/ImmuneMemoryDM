using Distributions
using SpecialFunctions
using JLD2
using FileIO
using Distributed
using LaTeXStrings
using Glob
using QuadGK


#This has the exponent (alpha d)^phi

function norm_gen_gaus(ϕ,α)
    return (2 /α * gamma(1+ 1/ϕ))
end


function integral_KB_del_ph(δ,ϕ,α)
    integral, err = quadgk(x -> exp(-α^(ϕ) *(abs(x))^ϕ)* (-α^(ϕ)* (abs((x)^ϕ) -abs((x-δ)^ϕ) )), - Inf, Inf, rtol=1e-3)
    return integral/norm_gen_gaus(ϕ,α)
end



function takememory_phi(β,v,α,ϕ)
    prob =probmemory_phi(β,v,α,ϕ)
    wahrscheinlikeit = rand()
    result=false
    if wahrscheinlikeit < prob
        result=true
    end
    result
end





function probmemory_phi(β,v,α,ϕ)
    prob = 1- exp(-β *energiememory_phi(v,α,ϕ) )
    prob
end

function energiememory_phi(v,α,ϕ)
    energie =exp(-(α^(ϕ)*(abs(v))^ϕ))/norm_gen_gaus(ϕ,α)
end



function freeenergy_real_phi(β,v,α,ϕ)
    energy = log( exp(β *energiememory_phi(v,α,ϕ)) -1)/β
    energy
end

function dissipation_phi(β,v,α,ϕ)
    diss = integral_KB_del_ph(v,ϕ,α)/β
    diss
end

function netUtility_phi(β,v,α,ϕ)
    unet=freeenergy_real_phi(β,v,α,ϕ)- dissipation_phi(β,v,α,ϕ)
    unet
end



function newv_drift(δ,n)
    vnew = sqrt(n)δ+randn()*δ*0.05*sqrt(n)
    #vnew = sqrt(n)δ+randn()*δ*0.1
    vnew
end



function onestep_dritf_real_phi(β,δ,oldv,α,n,ϕ)

    if (takememory_phi(β,oldv,α,ϕ))
        free=freeenergy_real_phi(β,oldv,α,ϕ)
        diss=dissipation_phi(β,oldv,α,ϕ)
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

function livesycle_drift_real_phi(β,δ,α,meanencounters,reps,ϕ)
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

             (netu,v,naiveuse,free,diss,n)=onestep_dritf_real_phi(β,δ,v,α,n,ϕ)
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








function findoptalphaoflive_drift_real_phi(β,δ,meanencounters,reps,αsteps,αmax,ϕ)

    start= 0.01
    α = range(start,stop=αmax,length=αsteps)

    data=livesycle_drift_real_phi.(β,δ,α,meanencounters,reps,ϕ)
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



function scandeltaliveGeqBparallel_drift_real_phi(β,meanencounters,reps,αsteps,αmax,deltasteps,ϕ)

    δ = range(0.005,stop=1.5/αmax,length= deltasteps)
    memstepsatmax =zeros(deltasteps)
    utilityopt =zeros(deltasteps)
    alphaopt =zeros(deltasteps)
    posnmax =zeros(deltasteps)
    freemax =zeros(deltasteps)
    dissmax =zeros(deltasteps)
    Threads.@threads for i =1:deltasteps
        memstepsatmax[i],utilityopt[i],alphaopt[i],posnmax[i],freemax[i],dissmax[i]=findoptalphaoflive_drift_real_phi(β,δ[i],meanencounters,reps,αsteps,αmax,ϕ)
    end
    return memstepsatmax,utilityopt,alphaopt,posnmax,δ,freemax,dissmax

end


function scanbetaeqGdeltaliveparallel_drift_real_phi(meanencounters,reps,αsteps,αmax,deltasteps,betasteps,ϕ)

    β = range(0.1,length=betasteps,stop=10.0)

    data = scandeltaliveGeqBparallel_drift_real_phi.(β,meanencounters,reps,αsteps,αmax,deltasteps,ϕ)

    memstepsatmax= zeros(betasteps,deltasteps)
    utilityopt= zeros(betasteps,deltasteps)
    alphaopt= zeros(betasteps,deltasteps)
    posnmax= zeros(betasteps,deltasteps)
    freemax =zeros(betasteps,deltasteps)
    dissmax =zeros(betasteps,deltasteps)
    for ii in 1:betasteps
        memstepsatmax[ii,:]=data[ii][1]
        utilityopt[ii,:]=data[ii][2]
        alphaopt[ii,:] = data[ii][3]
        posnmax[ii,:]=data[ii][4]
        freemax[ii,:]=data[ii][6]
        dissmax[ii,:]=data[ii][7]

    end
    δ = data[1][5]

    introduction="This gives the data meory_used_at_the_max  utility_atmax  otimal_alpha posible_maximum_memoryuses delta and beta and free_at_the_max and diss_at_the_max the data is Gamma first and delta second index"
 return memstepsatmax,utilityopt,alphaopt,posnmax,δ,β,freemax,dissmax,introduction

end












####Data Production



function simulation_DM_cluster_drift_real_oneB_phi(meanencounters,reps,asteps,bsteps,dsteps,amax,thisb,ϕ)

    βas = range(0.1,length=bsteps,stop=10.0)
    β=βas[thisb]
    parameters = (meanencounters,reps,asteps,bsteps,dsteps,amax,β,ϕ)


    filename = "Data/Gen_Bell/Singel_beta/meanencounters_$(meanencounters)/reps_$(reps)/Phi_$(ϕ)/DATA_Asteps_$(asteps)-Bsteps_$(bsteps)-Dsteps_$(dsteps)-amax_$(amax)_real-thisb_$(thisb).jld2"
    println("I want to do  $(filename)")
    if isfile(filename)
            println("Data $(filename) allready exists")
    else
        if !isdir("Data")
            mkdir("Data")
        end
        if !isdir("Data/Gen_Bell")
            mkdir("Data/Gen_Bell")
        end
        if !isdir("Data/Gen_Bell/Singel_beta")
            mkdir("Data/Gen_Bell/Singel_beta")
        end

        if !isdir("Data/Gen_Bell/Singel_beta/meanencounters_$(meanencounters)")
            mkdir("Data/Gen_Bell/Singel_beta/meanencounters_$(meanencounters)")
        end
        if !isdir("Data/Gen_Bell/Singel_beta/meanencounters_$(meanencounters)/reps_$(reps)")
            mkdir("Data/Gen_Bell/Singel_beta/meanencounters_$(meanencounters)/reps_$(reps)")
        end
        if !isdir("Data/Gen_Bell/Singel_beta/meanencounters_$(meanencounters)/reps_$(reps)/Phi_$(ϕ)")
            mkdir("Data/Gen_Bell/Singel_beta/meanencounters_$(meanencounters)/reps_$(reps)/Phi_$(ϕ)")
        end

        memstepsatmax,utilityopt,alphaopt,posnmax,δ,freemax,dissmax = scandeltaliveGeqBparallel_drift_real_phi(β,meanencounters,reps,asteps,amax,dsteps,ϕ)

    #filename = "/Users/oskar/Dropbox/PHD/Simulations/b-cell/Simulation/src/julia/data/test.jld"
        save(filename,"meory_used_at_the_min",  memstepsatmax,"max_utility",utilityopt,"otimal_alpha",alphaopt,"posible_maximum_memoryuses",posnmax,"parameters",parameters,"beta",β,"delta",δ,"freemax",freemax,"dissmax",dissmax,"thisb",thisb)
        println("I'm saving the file $(filename)")
    end
end





function addup_simulation_DM_cluster_drift_real_oneB_phi(meanencounters,reps,asteps,bsteps,dsteps,amax,ϕ)

    filename = "Data/Gen_Bell/full/meanencounters_$(meanencounters)/reps_$(reps)/Phi_$(ϕ)/DATA_Asteps_$(asteps)-Bsteps_$(bsteps)-Dsteps_$(dsteps)-amax_$(amax)_real.jld2"
    deltas = range(0.0001,length=dsteps,stop=0.4)

    parameters = (meanencounters,reps,asteps,bsteps,dsteps,amax,ϕ)



    if !isdir("Data/Gen_Bell/full")
        mkdir("Data/Gen_Bell/full")
    end
        if !isdir("Data/Gen_Bell/full/meanencounters_$(meanencounters)")
            mkdir("Data/Gen_Bell/full/meanencounters_$(meanencounters)")
        end
        if !isdir("Data/Gen_Bell/full/meanencounters_$(meanencounters)/reps_$(reps)")
            mkdir("Data/Gen_Bell/full/meanencounters_$(meanencounters)/reps_$(reps)")
        end
        if !isdir("Data/Gen_Bell/full/meanencounters_$(meanencounters)/reps_$(reps)/Phi_$(ϕ)")
            mkdir("Data/Gen_Bell/full/meanencounters_$(meanencounters)/reps_$(reps)/Phi_$(ϕ)")
        end



        βas = range(0.1,length=bsteps,stop=10.0)


    direct= "Data/Gen_Bell/Singel_beta/meanencounters_$(meanencounters)/reps_$(reps)/Phi_$(ϕ)"
    filenames = "DATA_Asteps_$(asteps)-Bsteps_$(bsteps)-Dsteps_$(dsteps)-amax_$(amax)_real-thisb_*.jld2"

    files2 = get_files(direct,filenames)
    l=length(files2)
    if l<bsteps
        println("!!!!  Achtung !!!!!!!!!!")
        println("es gibt nicht alle files ")
        println("!!!!  Achtung !!!!!!!!!!")
        println(" ")
        println("I have $(l)  of $(bsteps)")
        println(" ")
        println("Please do again ")
        println(" meanencounters_$(meanencounters)/reps_$(reps)/Phi_$(ϕ)/DATA_Asteps_$(asteps)-Bsteps_$(bsteps)-Dsteps_$(dsteps)-amax_$(amax)")
        #return -1
    end


    memstepsatmax= zeros(bsteps,dsteps)
    utilityopt= zeros(bsteps,dsteps)
    alphaopt= zeros(bsteps,dsteps)
    posnmax= zeros(bsteps,dsteps)
    freemax =zeros(bsteps,dsteps)
    dissmax =zeros(bsteps,dsteps)
    for ii in 1:length(files2)
        thisb = load(files2[ii],"thisb")
        memstepsatmax[thisb,:]=load(files2[ii],"meory_used_at_the_min")
        utilityopt[thisb,:]=load(files2[ii],"max_utility")
        alphaopt[thisb,:] = load(files2[ii],"otimal_alpha")
        posnmax[thisb,:]=load(files2[ii],"posible_maximum_memoryuses")
        freemax[thisb,:]=load(files2[ii],"freemax")
        dissmax[thisb,:]=load(files2[ii],"dissmax")
    end

    δ = load(files2[1],"delta")





    #filename = "/Users/oskar/Dropbox/PHD/Simulations/b-cell/Simulation/src/julia/data/test.jld"
        save(filename,"meory_used_at_the_min",  memstepsatmax,"max_utility",utilityopt,"otimal_alpha",alphaopt,"posible_maximum_memoryuses",posnmax,"parameters",parameters,"beta",βas,"delta",δ,"freemax",freemax,"dissmax",dissmax,"phi",ϕ)
            println("I'm saving the file $(filename)")

end



function get_files(direct,filenames)
    filenames=glob(filenames,direct)
    return filenames
end
