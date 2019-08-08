
###################################################
######### SSA ALGORITHM WITH PATH-MI   ############

function SSA(S::System,    # Variable storing system parameters and initial condition of species
            times::Vector{Float64};    # Vector of timepoints where to compute Path-MI . First entry is initial-condition time
            output_path::Bool=true, maxiters::Int64=200000,   # If and How much to pre-allocate to store the entire SSA path
            NstepMin::Int64=20,maxstep::Float64=0.01)         # Adaptive Euler integration settings
# Usage example:
# pmi, Xt = PathMI.SSA(S,times)
# pmi = PathMI.SSA(S,times, output_path=false)

# Initialize simulation variables
N,NReac=size(S.V)     # get number of chemical species and of reaction channels
@assert length(S.x0)==N "Initial condition doesn't match the dimensionality of the problem"
T=length(times)
time,tmax=times[1],times[end]
x=copy(S.x0)         # current value of state vector
output_path ? Xt=Path(time,S.x0,maxiters) : nothing  # initialize path object only if required by output_path
MomBgivenA=float([x[2];x[2].^2])  # zero-variance initial condition
MomAgivenB=float([x[1];x[1].^2])  # zero-variance initial condition
dMom=zeros(2)                     # dummy variable to work efficiently with conditional moments
MutualInfo=zeros(T)          # to store the pmi sample at the time-points defined in "times"
MIt=0.                       # current value of the pmi sample
prop=zeros(NReac)            # vector to compute SSA propensities
itercount=2
tcount=2
flag=true
while flag
    # Evaluate propensity functions
    compute_propensities(S.C,prop,x)
    PropSum=sum(prop)
    # Compute time of next reaction
    PropSum>0. ? WaitingTime=-1/PropSum*log(1-rand()) : WaitingTime=Inf
    timeold=time
    time+=WaitingTime
    # For all timepoints falling before the next event, evolve latent moments and output pmi at those times
    while times[tcount] < time
        deltat=times[tcount]-timeold
        ### until the new event, the SSA propensities are constant. Update directly their contribution to pmi
        MIt-=deltat*(prop[4]+prop[5]+prop[6])
        ### now evolve conditional moments and marginalized propensities
        MIt=EvolveAll(deltat,MIt,MomBgivenA,MomAgivenB,dMom,S.C,x,NstepMin,maxstep)
        MutualInfo[tcount]=MIt
        timeold=times[tcount]
        tcount+=1
        tcount>T ? (flag=false;time=tmax;break) : nothing
    end
    ### timeold is telling until where the pmi and the conditional moments have been propagated
    WaitingTime=time-timeold # note that it doesn't change if nothing was stored
    # If the next event happens before final time "tmax", update the system
    if time < tmax
        # Update moments and pmi until time
        MIt-=WaitingTime*(prop[4]+prop[5]+prop[6])
        MIt=EvolveAll(WaitingTime,MIt,MomBgivenA,MomAgivenB,dMom,S.C,x,NstepMin,maxstep)
        # Draw index of next reaction
        index=climbtower(rand()*PropSum,prop)
        # Update with jump terms if "index" is a coupling rection
        # PathMI has to be updated BEFORE the variables jump (cadlag!)
        if index==4
            MIt+=log(prop[4])-log(S.C[4]*x[1]*MomBgivenA[1])
            JumpMom(dMom,MomBgivenA)
        elseif index==5
            MIt+=log(prop[5])-log(S.C[5]*MomAgivenB[1])
            JumpMom(dMom,MomAgivenB)
        elseif index==6
            MIt+=log(prop[6])-log(S.C[6]*MomAgivenB[1]*x[2])
            JumpMom(dMom,MomAgivenB)
        end
        # Update state with chosen reaction
        for l=1:N
            x[l]+=S.V[l,index]
        end
        # Save path, if required
        if output_path
            if itercount==maxiters
                Xt.status=false
                break
            else
                save_path(itercount,Xt,index,time,x,MomAgivenB,MomBgivenA,MIt)
                itercount+=1
            end
        end
    end
end
if output_path
    Xt.status==false ? println("Overflow Error: Increase allocation size 'maxiters'! Break time: $time") : nothing
    save_path(itercount,Xt,0,tmax,x,MomAgivenB,MomBgivenA,MIt)
    cut_path(Xt,itercount)    # remove extra allocated entries in Path object
    return MutualInfo, Xt     # returns both pmi at timepoints="times" and Path of the system if "output_path"=true
else
    return MutualInfo         # returns only pmi at timepoints="times" if "output_path"=false
end
end  # END function



###################################################
############## ESTIMATION OF PATH-MI ##############

function estimate(S::System,   # Variable storing system parameters and initial condition of species
                  times::Vector{Float64},  # Vector of timepoints where to compute Path-MI . First entry is initial-condition time
                  Nsamples::Int64)  # How many Monte Carlo samples to use to estimate Path-MI
# Usage example:
# PMI, PMIstd, PMIrate, PMIratestd = PathPathMI.estimate(S,times,Nsamples)

InfoSum=zeros(length(times))
InfoSum2=copy(InfoSum)
println("You are computing the Path Mutual Information for the following network\n")
showSettings(S,times,Nsamples)
for i=1:Nsamples
    Info=SSA(S,times,output_path=false)
    InfoSum.+=Info
    InfoSum2.+=Info.^2
    #mod(i,100)==0 ? println("Iteration $i / $Nsamples") : nothing
end
InfoAv=InfoSum/Nsamples
InfoStdAv=sqrt.(InfoSum2/(Nsamples-1)-(InfoSum.^2)/(Nsamples*(Nsamples-1)))/sqrt(Nsamples)
InfoRateAv=InfoAv./times; InfoRateAv[1]=0. # correct 0./0.=NaN
InfoRateStdAv=InfoStdAv./times; InfoRateStdAv[1]=0. # correct 0./0.=NaN
return InfoAv, InfoStdAv, InfoRateAv, InfoRateStdAv
end


###################################################
####### FUNCTIONS GOVERNING TIME EVOLUTION ########

# COMPUTE PROPENSITY FUNCTIONS
function compute_propensities(C::Vector{Float64},a::Vector{Float64},x::Vector{Int64})
    # mass-action propensities for the system defined in create_network()
    a[1]=C[1]
    a[2]=C[2]*x[1]
    a[3]=C[3]*x[1]
    a[4]=C[4]*x[1]*x[2]
    a[5]=C[5]*x[1]
    a[6]=C[6]*x[1]*x[2]
    a[7]=C[7]*x[2]
return
end

# SYSTEM DEFINING EVOLUTION OF PATH-MI AND CONDITIONAL MOMENTS
function EvolveAll(tspan::Float64,MIt::Float64,MomBgivenA::Vector{Float64},MomAgivenB::Vector{Float64},dMom::Vector{Float64},C::Vector{Float64},x::Vector{Int64},NstepMin::Int64,maxstep::Float64)
# set adaptive Euler integration parameters
step,Nstep=integrationstep(tspan,NstepMin,maxstep)
for l=1:Nstep
    # first, compute marginalized propensities and update the Path-MI
    lambda4=C[4]*x[1]*MomBgivenA[1]
    MIt+=step*check_positive(lambda4)
    lambda5=C[5]*MomAgivenB[1]
    MIt+=step*check_positive(lambda5)
    lambda6=C[6]*MomAgivenB[1]*x[2]
    MIt+=step*check_positive(lambda6)
    # evolve the conditional moments, for both subnetwork selections
    EvolveBgivenA(step,dMom,MomBgivenA,C,x[1])
    EvolveAgivenB(step,dMom,MomAgivenB,C,x[2])
end
return MIt
end

# Evolves conditional moments of B when subnetwork A is selected
function EvolveBgivenA(step::Float64,dM::Vector{Float64},M::Vector{Float64},C::Vector{Float64},xa::Int64)
    bB=C[5]*xa
    bBA=C[6]*xa
    dM[1]=bB+(bBA-C[7])*M[1]-C[4]*xa*(M[2]-M[1]^2)
    dM[2]=bB+(2*bB+bBA+C[7])*M[1]+2*(bBA-C[7])*M[2]-2*C[4]*xa*M[2]/M[1]*(M[2]-M[1]^2) # Gamma closure
    M[1]+=step*dM[1]
    M[2]+=step*dM[2]
return
end

# Evolves conditional moments of A when subnetwork B is selected
function EvolveAgivenB(step::Float64,dM::Vector{Float64},M::Vector{Float64},C::Vector{Float64},xb::Int64)
    deathA=C[3]+C[4]*xb
    couplingB=C[5]+C[6]*xb
    dM[1]=C[1]+(C[2]-deathA)*M[1]-couplingB*(M[2]-M[1]^2)
    dM[2]=C[1]+(2*C[1]+C[2]+deathA)*M[1]+2*(C[2]-deathA)*M[2]-2*couplingB*M[2]/M[1]*(M[2]-M[1]^2) # Gamma closure
    M[1]+=step*dM[1]
    M[2]+=step*dM[2]
return
end

# Jump update the conditional moments when coupling reaction fires
function JumpMom(Mnew::Vector{Float64},M::Vector{Float64})
    Mnew[1]=M[2]/M[1]
    Mnew[2]=2*(M[2]/M[1])^2-M[2] ##M[2]+2*M[2]*(M[2]/M[1]^2-1) # Gamma closure
    M[1]=Mnew[1]
    M[2]=Mnew[2]
return
end

# Set adaptive Euler integration step
function integrationstep(wait::Float64,NstepMin::Int64,maxstep::Float64)
step=wait/NstepMin
if step>maxstep
    Nstep=floor(Int64,wait/maxstep)
    step=wait/Nstep
    return step,Nstep
else
    return step,NstepMin
end
end

# function to avoid to count NaN
check_positive(x::Float64)= x>0. ? (return x) : (return 0.)

# Perform Tower sampling
function climbtower(rr::Float64,vect::Vector{Float64})
i=1
cumul=vect[1]
while rr>cumul
   i+=1
   cumul+=vect[i]
end
return i
end


function showSettings(S::System, times::Vector{Float64}, Nsamples::Int64)
show(S)
println("from initial time $(times[1]) to final time $(times[end])")
println("and number of samples $Nsamples \n")
SSA(S,times,output_path=false) # discard first time
tc=@elapsed for i=1:10 SSA(S,times,output_path=false) end
println("On your machine, this computation is going to take approximately")
println("$(ceil(Int64,Nsamples*tc/10)) SECONDS\n\n")
return
end
