
font1 = Dict("family"=>"sans-serif","color"=>"black","weight"=>"normal","size"=>14)

# CASE STUDY 1: feed-forward catalytic mRNA->Protein network
function figure1a(;paper=false)
# Usage examples:
# timepoints, Amean, Bmean = launch.figure1a(paper=true)
# timepoints, Amean, Bmean = launch.figure1a()
    S=PathMI.initialize1()
    times=collect(0:1.:200)
    cc=["orange","green"]
    paper ? figure("figure1a",figsize=(4.2,4)) : figure("figure1a")
    paper ? Random.seed!(111) : nothing
    for i=1:10
        I,Xt=PathMI.SSA(S,times)
        plot(Xt.t,Xt.x[1,:],drawstyle="steps-post",color=cc[1],alpha=0.2,lw=0.5)
        plot(Xt.t,Xt.x[2,:],drawstyle="steps-post",color=cc[2],alpha=0.2,lw=0.5)
    end
    Amean=S.C[1]/S.C[3].+(S.x0[1].-S.C[1]/S.C[3])*exp.(-S.C[3]*times)
    Bmean=S.C[5]*Amean/S.C[7].+(S.x0[2].-S.C[5]*Amean/S.C[7]).*exp.(-S.C[7]*times)
    plot(times,Bmean,color=cc[2],label="B")
    plot(times,Amean,color=cc[1],label="A")  # already initialized at equilibrium
    legend()
    xlabel("time",fontdict=font1)
    ylabel("copy number",fontdict=font1)
    paper ? subplots_adjust(top=0.995,bottom=0.115,left=0.17,right=0.995,hspace=0.2,wspace=0.2) : nothing
return times,Amean,Bmean
end

function figure1b(;paper=false)
# Usage examples:
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure1b(paper=true)
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure1b()
    S=PathMI.initialize1()
    times=collect(0:1.:200)
    Nsamples=10000
    paper ? Random.seed!(111) : nothing
    # Compute Path Mutual Information
    PMI,PMIstd, PMIrate, PMIratestd=PathMI.estimate(S,times,Nsamples)
    # Plot Path Mutual Information
    paper ? figure("figure1b",figsize=(4,4)) : figure("figure1b")
    plot(times,PMI,color="xkcd:royal blue")
    xlabel("time",fontdict=font1)
    ylabel("path mutual information [nats]",fontdict=font1)
    # Plot single realizations
    for i=1:20
        I,Xt=PathMI.SSA(S,times)
        figure("figure1b")
        plot(Xt.t,Xt.pmi,color="xkcd:royal blue",alpha=0.2,lw=0.5)
    end
    paper ? (figure("figure1b"); subplots_adjust(top=0.995,bottom=0.115,left=0.13,right=0.995,hspace=0.2,wspace=0.2) ) : nothing
return times, S, PMI, PMIstd, PMIrate, PMIratestd
end


# CASE STUDY 2: Transcription burst mRNA->Protein
function figure1c(;paper=false)
    # Usage examples:
    # timepoints, Amean, Bmean = launch.figure1c(paper=true)
    # timepoints, Amean, Bmean = launch.figure1c()
    S=PathMI.initialize2()
    times=collect(0:1.:200)
    cc=["orange","green"]
    paper ? figure("figure1c",figsize=(4.2,4)) : figure("figure1c")
    paper ? Random.seed!(222) : nothing
    for i=1:10
        I,Xt=PathMI.SSA(S,times)
        plot(Xt.t,Xt.x[1,:],drawstyle="steps-post",color=cc[1],alpha=0.2,lw=0.5)
        plot(Xt.t,Xt.x[2,:],drawstyle="steps-post",color=cc[2],alpha=0.2,lw=0.5)
    end
    Amean=S.C[1]/S.C[3].+(S.x0[1].-S.C[1]/S.C[3])*exp.(-S.C[3]*times)
    S.C[3]==S.C[7] ? ct=times.*exp.(-S.C[7]*times) : ct=(exp.(-S.C[3]*times).-exp.(-S.C[7]*times))/(S.C[7]-S.C[3])
    Bmean=(S.C[5]*S.C[1])/(S.C[7]*S.C[3]).+(S.x0[2]-(S.C[5]*S.C[1])/(S.C[7]*S.C[3]))*exp.(-S.C[7]*times).+S.C[5]*(S.x0[1]-S.C[1]/S.C[3]).*ct
    plot(times,Bmean,color=cc[2],label="B")
    plot(times,Amean,color=cc[1],label="A")  # already initialized at equilibrium
    legend()
    xlabel("time",fontdict=font1)
    ylabel("copy number",fontdict=font1)
    paper ? subplots_adjust(top=0.995,bottom=0.115,left=0.17,right=0.995,hspace=0.2,wspace=0.2) : nothing
return times,Amean,Bmean
end

function figure1d(;paper=false)
# Usage examples:
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure1d(paper=true)
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure1d()
    S=PathMI.initialize2()
    times=collect(0:1.:200)
    Nsamples=10000
    paper ? Random.seed!(1234567) : nothing
    # Compute Path Mutual Information
    PMI,PMIstd, PMIrate, PMIratestd=PathMI.estimate(S,times,Nsamples)
    # Plot Path Mutual Information
    paper ? figure("figure1d",figsize=(4.,4)) : figure("figure1d")
    plot(times,PMI,color="xkcd:royal blue")
    xlabel("time",size=12)
    ylabel("path mutual information [nats]",fontdict=font1)
    # Plot single realizations
    for i=1:20
        I,Xt=PathMI.SSA(S,times)
        # in Xt only event-times are stored. For better visualization, join also the portion of PMI in I after species get extincted
        index_A0=findall(Xt.x[1,:].==0)[1]
        indexes_to_join=findall(times.>Xt.t[index_A0])
        tt=[Xt.t[1:index_A0];times[indexes_to_join]]
        Info=[Xt.pmi[1:index_A0];I[indexes_to_join]]
        figure("figure1d")
        plot(tt,Info,color="xkcd:royal blue",alpha=0.2,lw=0.5)
    end
    paper ? (figure("figure1d"); subplots_adjust(top=0.995,bottom=0.115,left=0.13,right=0.995,hspace=0.2,wspace=0.2) ) : nothing
return times, S, PMI, PMIstd, PMIrate, PMIratestd
end


# CASE STUDY 3: Lotka-Volterra model
function figure2a(;paper=false)
# Usage examples:
# timepoints, pmi = launch.figure2a(paper=true)
# timepoints, pmi = launch.figure2a()
    S=PathMI.initialize3()
    times=collect(0:1.:300)
    figure("figure2a",figsize=(7.2,4))
    cc=["orange","green"]
    paper ? Random.seed!(790966) : nothing
    # Compute one realization
    I,Xt=PathMI.SSA(S,times,maxiters=1000000)
    # Preparing to plot
    tt=Xt.t
    Info=Xt.pmi
    indexes_A0=findall(Xt.x[1,:].==0)
    if length(indexes_A0)>0
        index_A0=indexes_A0[1]
        indexes_to_join=findall(times.>Xt.t[index_A0])
        tt=[Xt.t[1:index_A0];times[indexes_to_join]]
        Info=[Xt.pmi[1:index_A0];I[indexes_to_join]]
    end
    # Plot species trajectories
    plot(Xt.t,Xt.x[2,:],color=cc[2],alpha=0.8,lw=1.,label="B")
    plot(Xt.t,Xt.x[1,:],color=cc[1],alpha=0.8,lw=1.,label="A")
    xlabel("time",fontdict=font1)
    ylabel("copy number",fontdict=font1)
    legend(loc="upper left")
    ax=gca()
    ax2 = ax.twinx() # Create another axis on top of the current axis
    # Plot PMI sample
    plot(tt,Info,color="xkcd:royal blue",alpha=0.8,lw=1.)
    paper ? ylim(nothing,150) : nothing
    setp(ax2.get_yticklabels(),color="xkcd:royal blue") # Y Axis font formatting
    ylabel("path-MI sample",size=14, color="xkcd:royal blue")
return tt,Info
end

function figure2b(;paper=false)
# Usage examples:
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure2b(paper=true)
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure2b()
    S=PathMI.initialize3()
    times=collect(0:1.:1000)
    Nsamples=1000
    paper ? Random.seed!(1234) : nothing
    # Compute Path Mutual Information
    PMI,PMIstd, PMIrate, PMIratestd=PathMI.estimate(S,times,Nsamples)
    # Plot Path Mutual Information
    figure("figure2b",figsize=(7,4))
    plot(times,PMI,color="xkcd:royal blue")
    xlabel("time",fontdict=font1)
    ylabel("Path mutual information [nats]",fontdict=font1)
    # Plot single realizations
    for i=1:20
        I,Xt=PathMI.SSA(S,times,maxiters=2000000)
        # in Xt only event-times are stored. For better visualization, join also the portion of PMI in I after species get extincted
        tt=Xt.t
        Info=Xt.pmi
        indexes_A0=findall(Xt.x[1,:].==0)
        if length(indexes_A0)>0
            index_A0=indexes_A0[1]
            indexes_to_join=findall(times.>Xt.t[index_A0])
            tt=[Xt.t[1:index_A0];times[indexes_to_join]]
            Info=[Xt.pmi[1:index_A0];I[indexes_to_join]]
        end
        figure("figure2b")
        plot(tt,Info,color="xkcd:royal blue",alpha=0.2,lw=0.5)
        paper ? subplots_adjust(left=0.11,right=0.91) : nothing
    end
return times, S, PMI, PMIstd, PMIrate, PMIratestd
end
