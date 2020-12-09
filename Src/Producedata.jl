using Distributions
using SpecialFunctions
using JLD2
using FileIO
using LaTeXStrings
using Glob


    include("simulations-bcel-DM-diss.jl")


function simulation_oneB(meanencounters,reps,asteps,bsteps,dsteps,amax,thisb)

    βas = range(0.1,length=bsteps,stop=10.0)
    β=βas[thisb]
    parameters = (meanencounters,reps,asteps,bsteps,dsteps,amax,β)


    filename = "Data/Drift_data/SingelBs/meanencounters_$(meanencounters)/reps_$(reps)/DATA_Asteps_$(asteps)-Bsteps_$(bsteps)-Dsteps_$(dsteps)-amax_$(amax)_real-thisb_$(thisb).jld2"
    println("I want to do  $(filename)")
    if isfile(filename)
            println("Data $(filename) allready exists")
    else
        if !isdir("Data")
            mkdir("Data")
        end
        if !isdir("Data/Drift_data")
            mkdir("Data/Drift_data")
        end
        if !isdir("Data/Drift_data/SingelBs")
            mkdir("Data/Drift_data/SingelBs")
        end

        if !isdir("Data/Drift_data/SingelBs/meanencounters_$(meanencounters)")
            mkdir("Data/Drift_data/SingelBs/meanencounters_$(meanencounters)")
        end
        if !isdir("Data/Drift_data/SingelBs/meanencounters_$(meanencounters)/reps_$(reps)")
            mkdir("Data/Drift_data/SingelBs/meanencounters_$(meanencounters)/reps_$(reps)")
        end


        memstepsatmax,utilityopt,alphaopt,posnmax,δ,freemax,dissmax = find_optimal_alpha_for_range_of_delta(β,meanencounters,reps,asteps,amax,dsteps)

    #filename = "/Users/oskar/Dropbox/PHD/Simulations/b-cell/Simulation/src/julia/data/test.jld"
        save(filename,"meory_used_at_the_min",  memstepsatmax,"max_utility",utilityopt,"otimal_alpha",alphaopt,"posible_maximum_memoryuses",posnmax,"parameters",parameters,"beta",β,"delta",δ,"freemax",freemax,"dissmax",dissmax,"thisb",thisb)
        println("I'm saving the file $(filename)")
    end
end





function addup_oneB(meanencounters,reps,asteps,bsteps,dsteps,amax)

    filename = "Data/Drift_data/Fulldata/meanencounters_$(meanencounters)/reps_$(reps)/DATA_Asteps_$(asteps)-Bsteps_$(bsteps)-Dsteps_$(dsteps)-amax_$(amax)_real.jld2"
    println(filename)
    if !isdir("Data")
        mkdir("Data")
    end
    if !isdir("Data/Drift_data")
        mkdir("Data/Drift_data")
    end
    if !isdir("Data/Drift_data/Fulldata")
        mkdir("Data/Drift_data/Fulldata")
    end

        if !isdir("Data/Drift_data/Fulldata/meanencounters_$(meanencounters)")
            mkdir("Data/Drift_data/Fulldata/meanencounters_$(meanencounters)")
        end
        if !isdir("Data/Drift_data/Fulldata/meanencounters_$(meanencounters)/reps_$(reps)")
            mkdir("Data/Drift_data/Fulldata/meanencounters_$(meanencounters)/reps_$(reps)")
        end



    βas = range(0.1,length=bsteps,stop=10.0)
    parameters = (meanencounters,reps,asteps,bsteps,dsteps,amax,βas)

    direct= "Data/Drift_data/SingelBs/meanencounters_$(meanencounters)/reps_$(reps)/"
    filenames = "DATA_Asteps_$(asteps)-Bsteps_$(bsteps)-Dsteps_$(dsteps)-amax_$(amax)_real-thisb_*.jld2"
    files2 = get_files(direct,filenames)

    l=length(files2)
    println(l)
    if l<bsteps
        println("!!!!  Achtung !!!!!!!!!!")
        println("es gibt nicht alle files ")
        println("!!!!  Achtung !!!!!!!!!!")
        println(" ")
        println("I have $(l)  of $(bsteps)")
        println(" ")
        println("Please do again ")
        println(" meanencounters_$(meanencounters)/reps_$(reps)/DATA_Asteps_$(asteps)-Bsteps_$(bsteps)-Dsteps_$(dsteps)-amax_$(amax)")
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
        println("done")



println(filename)
        save(filename,"meory_used_at_the_min",  memstepsatmax,"max_utility",utilityopt,"otimal_alpha",alphaopt,"posible_maximum_memoryuses",posnmax,"parameters",parameters,"beta",βas,"delta",δ,"freemax",freemax,"dissmax",dissmax)
        println("I'm saving the file $(filename)")

end




# produce the distribution of alphas



function get_files(direct,filenames)
    filenames=glob(filenames,direct)
    return filenames
end
