module launch

using Random
using PyPlot

include("PathMI.jl")

font1 = Dict("family"=>"sans-serif","color"=>"black","weight"=>"normal","size"=>14)

function figure1b(;paper=false)
# feed-forward catalytic mRNA->Protein network
# feed-forward catalytic mRNA->Protein network
# Usage examples:
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure1b(paper=true)
### to get the identical plot of the paper
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure1b()
### to run new simulations
S=PathMI.initialize1()
times=collect(0:1.:200)
Nsamples=10000
paper ? Random.seed!(111) : nothing
PMI,PMIstd, PMIrate, PMIratestd=PathMI.estimate(S,times,Nsamples)
paper ? figure("figure1b",figsize=(4,4)) : figure("figure1b")
plot(times,PMI,color="xkcd:royal blue")
xlabel("time",fontdict=font1)
ylabel("path mutual information [nats]",fontdict=font1)
#paper ? figure("MIrate1",figsize=(4.2,4)) : figure("MIrate1")
#plot(times,PMIrate,color="xkcd:royal blue")
#xlabel("time",fontdict=font1)
#ylabel("path mutual information rate",fontdict=font1)
for i=1:20
    I,Xt=PathMI.SSA(S,times)
    plotpoints=plotindices(Xt.t)
    figure("figure1b")
    plot(Xt.t[plotpoints],Xt.pmi[plotpoints],color="xkcd:royal blue",alpha=0.2,lw=0.5)
    #figure("MIrate1")
    #plot(Xt.t[plotpoints],[0.;Xt.pmi[plotpoints[2:end]]]./[0.;Xt.t[plotpoints[2:end]]],color="xkcd:royal blue",alpha=0.1,lw=0.5)
end
if paper
    figure("figure1b")
    subplots_adjust(top=0.995,bottom=0.115,left=0.13,right=0.995,hspace=0.2,wspace=0.2)
    #figure("MIrate1")
    #subplots_adjust(top=0.995,bottom=0.115,left=0.17,right=0.995,hspace=0.2,wspace=0.2)
    #ylim(-0.25,0.75)
end
return times, S, PMI, PMIstd, PMIrate, PMIratestd
end


function figure1a(;paper=false)
    # feed-forward catalytic mRNA->Protein network
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


function figure1d(;paper=false)
# Transcription burst mRNA->Protein
# Transcription burst mRNA->Protein
# Usage examples:
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure1d(paper=true)
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure1d()
S=PathMI.initialize2()
times=collect(0:1.:200)
Nsamples=10000
paper ? Random.seed!(1234567) : nothing
PMI,PMIstd, PMIrate, PMIratestd=PathMI.estimate(S,times,Nsamples)
paper ? figure("figure1d",figsize=(4.,4)) : figure("figure1d")
plot(times,PMI,color="xkcd:royal blue")
xlabel("time",size=12)
ylabel("path mutual information [nats]",fontdict=font1)
#paper ? figure("MIrate2",figsize=(4.2,4)) : figure("MIrate2")
#plot(times,PMIrate,color="xkcd:royal blue")
#xlabel("time",size=12)
#ylabel("path mutual information rate [nats]",fontdict=font1)
for i=1:20
    I,Xt=PathMI.SSA(S,times)
    A0=findall(Xt.x[1,:].==0)[1]
    join=findall(times.>Xt.t[A0])
    tt=Xt.t[1:A0]
    Info=Xt.pmi[1:A0]
    plotpoints=plotindices(tt)
    tt=[tt[plotpoints];times[join]]
    Info=[Info[plotpoints];I[join]]
    figure("figure1d")
    plot(tt,Info,color="xkcd:royal blue",alpha=0.2,lw=0.5)
    #figure("MIrate2")
    #plot(tt,[0.;Info[2:end]./tt[2:end]],color="xkcd:royal blue",alpha=0.1,lw=0.5)
end
if paper
    figure("figure1d")
    subplots_adjust(top=0.995,bottom=0.115,left=0.13,right=0.995,hspace=0.2,wspace=0.2)
    #figure("MIrate2")
    #subplots_adjust(top=0.995,bottom=0.115,left=0.17,right=0.995,hspace=0.2,wspace=0.2)
    #ylim(-0.125,0.55)
end
return times, S, PMI, PMIstd, PMIrate, PMIratestd
end


function figure1c(;paper=false)
    # Transcription burst mRNA->Protein
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
        plotpoints=plotindices(Xt.t)
        plot(Xt.t[plotpoints],Xt.x[1,plotpoints],drawstyle="steps-post",color=cc[1],alpha=0.2,lw=0.5)
        plot(Xt.t[plotpoints],Xt.x[2,plotpoints],drawstyle="steps-post",color=cc[2],alpha=0.2,lw=0.5)
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

function figure2b(;paper=false)
# Lotka-Volterra model
# Usage examples:
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure2b(paper=true)
# timepoints, System, PMI, PMIstd, PMIrate, PMIratestd = launch.figure2b()
S=PathMI.initialize3()
times=collect(0:1.:1000)
Nsamples=1000
paper ? Random.seed!(1234) : nothing
PMI,PMIstd, PMIrate, PMIratestd=PathMI.estimate(S,times,Nsamples)
figure("figure2b",figsize=(7,4))
plot(times,PMI,color="xkcd:royal blue")
xlabel("time",fontdict=font1)
ylabel("Path mutual information [nats]",fontdict=font1)
#figure("MIrate3",figsize=(7.2,4))
#plot(times,PMIrate,color="xkcd:royal blue")
#xlabel("time",fontdict=font1)
#ylabel("Path mutual information rate",fontdict=font1)
for i=1:20
    I,Xt=PathMI.SSA(S,times,maxiters=2000000)
    A0=findall(Xt.x[1,:].==0)
    if length(A0)>0
        A0=A0[1]
        join=findall(times.>Xt.t[A0])
        tt=Xt.t[1:A0]
        Info=Xt.pmi[1:A0]
        plotpoints=plotindices(tt)
        tt=[tt[plotpoints];times[join]]
        Info=[Info[plotpoints];I[join]]
    else
        plotpoints=plotindices(Xt.t)
        tt=Xt.t[plotpoints]
        Info=Xt.pmi[plotpoints]
    end
    figure("figure2b")
    plot(tt,Info,color="xkcd:royal blue",alpha=0.2,lw=0.5)
    paper ? subplots_adjust(left=0.11,right=0.91) : nothing
    #figure("MIrate3")
    #plot(tt,[0.;Info[2:end]./tt[2:end]],color="xkcd:royal blue",alpha=0.1,lw=0.5)
    #paper ? (ylim(0.,1.45); yticks([0;0.25;0.5;0.75;1.;1.25])) : nothing
end
return times, S, PMI, PMIstd, PMIrate, PMIratestd
end


function figure2a(;paper=false)
    # Lotka-Volterra model
    # Usage examples:
    # timepoints, pmi = launch.figure2a(paper=true)
    # timepoints, pmi = launch.figure2a()
    S=PathMI.initialize3()
    times=collect(0:1.:300)
    figure("figure2a",figsize=(7.2,4))
    cc=["orange","green"]
    paper ? Random.seed!(790966) : nothing
    I,Xt=PathMI.SSA(S,times,maxiters=1000000)
    A0=findall(Xt.x[1,:].==0)
    if length(A0)>0
        A0=A0[1]
        join=findall(times.>Xt.t[A0])
        tt=Xt.t[1:A0]
        Info=Xt.pmi[1:A0]
        plotpoints=plotindices(tt)
        tt=[tt[plotpoints];times[join]]
        Info=[Info[plotpoints];I[join]]
    else
        plotpoints=plotindices(Xt.t)
        tt=Xt.t[plotpoints]
        Info=Xt.pmi[plotpoints]
    end
    plotpoints2=plotindices(Xt.t,NumberToPlot=5000)
    #subplot(2,1,1)
    plot(Xt.t[plotpoints2],Xt.x[2,plotpoints2],color=cc[2],alpha=0.8,lw=1.,label="B")
    plot(Xt.t[plotpoints2],Xt.x[1,plotpoints2],color=cc[1],alpha=0.8,lw=1.,label="A")
    xlabel("time",fontdict=font1)
    ylabel("copy number",fontdict=font1)
    legend(loc="upper left")
    ax=gca()
    ax2 = ax.twinx() # Create another axis on top of the current axis
    #subplot(2,1,2)
    plot(tt,Info,color="xkcd:royal blue",alpha=0.8,lw=1.)
    paper ? ylim(nothing,150) : nothing
    setp(ax2.get_yticklabels(),color="xkcd:royal blue") # Y Axis font formatting
    ylabel("path-MI sample",size=14, color="xkcd:royal blue")
return tt,Info
end

# function to use less points in plots not to make heavy vector figures
function plotindices(Xtt::Vector{Float64};NumberToPlot::Int64=1000)
    PP=length(Xtt)
    if PP>NumberToPlot
        skip=floor(Int64,PP/NumberToPlot)
        plotpoints=[1;collect(2:skip:PP)]
        plotpoints[end]==PP ? nothing : plotpoints=[plotpoints;PP]
        return plotpoints
    else
        return collect(1:length(Xtt))
    end
end

end  #END MODULE
