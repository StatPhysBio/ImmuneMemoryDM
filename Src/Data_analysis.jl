using Distributions
using SpecialFunctions
using JLD2
using FileIO
using Distributed










#### Analysis of double max


function build_double_opt_for_files(files,κ)
    num=length(files)
    δ= load(files[1],"delta")
    dsteps=length(δ)
    meanencouters_of_file=zeros(num)
    doubleopt=zeros(num,dsteps)
    aopt=zeros(num,dsteps)
    bopt=zeros(num,dsteps)
    dissopt=zeros(num,dsteps)
    usesopt=zeros(num,dsteps)
    for jj in 1:num
        δ= load(files[jj],"delta")
        β = load(files[jj],"beta")
        meanencounters2,repsOFproduction,asteps,bsteps,dsteps,amax= load(files[jj],"parameters")
        utilityopt1load = load(files[jj],"max_utility")
        utilityopt= load(files[jj],"max_utility")
        alphaopt= load(files[jj],"otimal_alpha")
        diss= load(files[jj],"dissmax")
        uses= load(files[jj],"meory_used_at_the_min")

        for ii in 1:dsteps
            for kk in 1:bsteps
                utilityopt[kk,ii] = utilityopt1load[kk,ii] -meanencounters2 *κ*β[kk]
            end
        end
        doubleopt[jj,:],aopt[jj,:],bopt[jj,:],dissopt[jj,:],usesopt[jj,:]=give_double_max_from_max_a_all(utilityopt./meanencounters2,alphaopt,β,diss./meanencounters2,uses)

        meanencouters_of_file[jj]=meanencounters2
    end
    return meanencouters_of_file,doubleopt,aopt,bopt,dissopt,usesopt,δ
end

function build_double_opt_for_files_amax(files,κ)
    num=length(files)
    δ= load(files[1],"delta")
    dsteps=length(δ)
    meanencouters_of_file=zeros(num)
    amaxes=zeros(num)
    doubleopt=zeros(num,dsteps)
    aopt=zeros(num,dsteps)
    bopt=zeros(num,dsteps)
    dissopt=zeros(num,dsteps)
    usesopt=zeros(num,dsteps)
    for jj in 1:num
        δ= load(files[jj],"delta")
        β = load(files[jj],"beta")
        meanencounters2,repsOFproduction,asteps,bsteps,dsteps,amax= load(files[jj],"parameters")
        utilityopt1load = load(files[jj],"max_utility")
        utilityopt= load(files[jj],"max_utility")
        alphaopt= load(files[jj],"otimal_alpha")
        diss= load(files[jj],"dissmax")
        uses= load(files[jj],"meory_used_at_the_min")

        for ii in 1:dsteps
            for kk in 1:bsteps
                utilityopt[kk,ii] = utilityopt1load[kk,ii] -meanencounters2 *κ*β[kk]
            end
        end
        doubleopt[jj,:],aopt[jj,:],bopt[jj,:],dissopt[jj,:],usesopt[jj,:]=give_double_max_from_max_a_all(utilityopt./meanencounters2,alphaopt,β,diss./meanencounters2,uses)
        amaxes[jj]=amax
        meanencouters_of_file[jj]=meanencounters2
    end
    return meanencouters_of_file,doubleopt,aopt,bopt,dissopt,usesopt,δ,amaxes
end





function give_double_max_from_max_a_all(utilityopt,alphaopt,β,diss,uses)
    ds=length(utilityopt[1,:])
    doubleopt=zeros(ds)
    aopt=zeros(ds)
    bopt=zeros(ds)
    dissopt=zeros(ds)
    usesopt=zeros(ds)
    for ii in 1:ds
        argbopt=argmax(utilityopt[:,ii])
        doubleopt[ii]=utilityopt[argbopt,ii]
        dissopt[ii]=diss[argbopt,ii]
        usesopt[ii]=uses[argbopt,ii]
        aopt[ii]=alphaopt[argbopt,ii]
        bopt[ii]=β[argbopt]
    end
    return     doubleopt,aopt,bopt,dissopt,usesopt
end






####File Analysis

function analysis_of_opt_for_kappa(filename,k_steps)
    κ=range(0,0.2,length=k_steps)
    utilityopt1load = load(filename,"max_utility")
    alphaopt2load= load(filename,"otimal_alpha")
    dissmax2load= load(filename,"dissmax")
    usesload= load(filename,"meory_used_at_the_min")
    δ= load(filename,"delta")
    β = load(filename,"beta")
    freeload=utilityopt1load.+dissmax2load
    meanencounters2,reps,asteps,bsteps,dsteps,amax= load(filename,"parameters")
    doubleopt=zeros(k_steps,dsteps)
    aopt=zeros(k_steps,dsteps)
    bopt=zeros(k_steps,dsteps)
    dissopt=zeros(k_steps,dsteps)
    usesopt=zeros(k_steps,dsteps)
    for kk in 1:k_steps
        utilityopt1 = zeros(bsteps,dsteps)
       for ii in 1:dsteps
            for jj in 1:bsteps
                 utilityopt1[jj,ii] = utilityopt1load[jj,ii] -meanencounters2 *κ[kk]*β[jj]
            end
        end
        doubleopt[kk,:],aopt[kk,:],bopt[kk,:],dissopt[kk,:],usesopt[kk,:]=give_double_max_from_max_a_all(utilityopt1,alphaopt2load,β,dissmax2load./meanencounters2,usesload)

    end
    return doubleopt,aopt,bopt,dissopt,usesopt, δ
end


function put_data_in_region_and_order(x,y,xlim,ylim)
    points=length(x)
    newpoints=0
    (xdown,xup)=xlim
    (ydown,yup)=ylim
    for ii in 1:points
        xi=x[ii]
        yi=y[ii]
        if (xi<xup && xi>xdown && yi<yup && yi>ydown)
            newpoints +=1
        end
    end
    xregion=zeros(newpoints)
    yregion=zeros(newpoints)
    jj=1
    for ii in 1:points
        xi=x[ii]
        yi=y[ii]
        if (xi<xup && xi>xdown && yi<yup && yi>ydown)
            xregion[jj]=xi *1
            yregion[jj]=yi*1
            jj +=1
        end
    end
     xorderd=zeros(newpoints)
    yorderd=zeros(newpoints)
    top=maximum(x)+10
    for jj in 1:newpoints
        argument= argmin(xregion)
        xorderd[jj]=xregion[argument]*1
        yorderd[jj]=yregion[argument]*1
        xregion[argument]=top*2
    end
    return xorderd,yorderd
end

function build_diff_meanencouters_for_files(files,meanencountersmax,reps,pos_delta,κ)
    num=length(files)
    utlity=zeros(num,meanencountersmax)
    dissopt=zeros(num,meanencountersmax)
    δ= load(files[1],"delta")
    meanencounters=range(1,length=meanencountersmax,step=1)
    meanencouters_of_file=zeros(num)
    uses=zeros(num,meanencountersmax)
    for jj in 1:num
        δ= load(files[jj],"delta")
        β = load(files[jj],"beta")
        meanencounters2,repsOFproduction,asteps,bsteps,dsteps,amax= load(files[jj],"parameters")
        utilityopt1load = load(files[jj],"max_utility")
        utilityopt= load(files[jj],"max_utility")
        alphaopt= load(files[jj],"otimal_alpha")
        for ii in 1:dsteps
            for kk in 1:bsteps
                utilityopt[kk,ii] = utilityopt1load[kk,ii] -meanencounters2 *κ*β[kk]
            end
        end
        doubleopt,aopt,bopt=give_double_max_from_max_a(utilityopt./meanencounters2,alphaopt,β)
        println("I found the opt its at $(bopt[pos_delta]), $(δ[pos_delta]) and $(aopt[pos_delta])")
        dummi0,dummi1,dummi2,dummi3,dummi4=test_dif_meanencouters(bopt[pos_delta],δ[pos_delta],aopt[pos_delta],meanencountersmax,reps)
        utlity[jj,:] = dummi0./dummi4 .- bopt[pos_delta] * κ
        dissopt[jj,:] = dummi3./dummi4
        uses[jj,:]=dummi1
        meanencouters_of_file[jj]=meanencounters2
    end
    return meanencouters_of_file,utlity,dissopt,meanencounters,uses
end


function test_dif_meanencouters(β,δ,α,meanencountersmax,reps)
    utilty=zeros(meanencountersmax)
    uses=zeros(meanencountersmax)
    free=zeros(meanencountersmax)
    diss=zeros(meanencountersmax)
    Threads.@threads for ii in 1:meanencountersmax
        utilty[ii],uses[ii],free[ii],diss[ii]=livesycle(β,δ,α,ii,reps)
    end
    meanencounters=range(1,length=meanencountersmax,step=1)
    return utilty,uses,free,diss,meanencounters
end


function builed_transistions_b_scaled(files,k_steps)
    κ=range(0,0.01/4,length=k_steps)
    κ=κ.*4/sqrt(π)
    transision_threshold=0.4 # persent of change which is a transition
     meanencouters_of_file,doubleopt1,aopt1,bopt1,dissopt1,usesopt1,δ= build_double_opt_for_files(files,0)
    dsteps=length(δ)
    encounter_steps=length(meanencouters_of_file)
    trans_points=zeros(encounter_steps,k_steps)
    for kk in 1:k_steps
        meanencouters_of_file,doubleopt,aopt,bopt,dissopt,usesopt,δ= build_double_opt_for_files(files,κ[kk])
        for ee in 1:encounter_steps
            ddiss=0.0
            transdiss=1
            for dd in 1:dsteps-1

                bdiss=(-bopt[ee,dd+1]+bopt[ee,dd])/bopt[ee,dd+1]
                if ((bdiss>transision_threshold) && bopt[ee,dd+1]>1)

                    transdiss=dd
                end
            end
            trans_points[ee,kk]=δ[transdiss]

        end
    end
    return meanencouters_of_file,δ,κ,trans_points
end


function give_double_max_from_max_a(utilityopt,alphaopt,β)
    ds=length(utilityopt[1,:])
    doubleopt=zeros(ds)
    aopt=zeros(ds)
    bopt=zeros(ds)
    for ii in 1:ds
        argbopt=argmax(utilityopt[:,ii])
        doubleopt[ii]=utilityopt[argbopt,ii]
        aopt[ii]=alphaopt[argbopt,ii]
        bopt[ii]=β[argbopt]
    end
    return     doubleopt,aopt,bopt
end
