
# OBJECT TO STORE SYSTEM PARAMETERS
mutable struct System
    C::Vector{Float64}    # Vector of the rate constants
    V::Matrix{Int64}      # Matrix of stochiometric coefficients
    x0::Vector{Int64}     # Initial condition
    function System(C::Vector{Float64},V::Matrix{Int64},x0::Vector{Int64})
        @assert length(x0)==size(V,1)    "Dimensions of initial condition and stoichiometric matrix do not match"
        @assert length(C)==size(V,2)     "Dimensions of reaction rates and stoichiometric matrix do not match"
    return new(C,V,x0)
    end
end

function create_network()
# IMPLEMENTED REACTIONS
#1# 0->A
#2# A->2A
#3# A->0
#4# A+B->B
#5# A->A+B
#6# A+B->A+2B
#7# B->0
# State-change Matrix :
   #1# #2# #3#  #4#  #5#   #6#    #7#
V=[1    1  -1   -1    0     0     0;  # A
   0    0   0    0    1     1    -1]  # B
C=zeros(size(V,2))          # Rate constants vector
x0=zeros(Int64,size(V,1))   # Initial condition
return System(C,V,x0)
end

function Base.show(io::IO, S::System)
S.C[1]>0. ? println("0 -> A     with rate  $(S.C[1])") : nothing
S.C[2]>0. ? println("A -> 2A     with rate  $(S.C[2])") : nothing
S.C[3]>0. ? println("A -> 0     with rate  $(S.C[3])") : nothing
S.C[4]>0. ? println("A + B -> B      with rate  $(S.C[4])") : nothing
S.C[5]>0. ? println("A -> A + B     with rate  $(S.C[5])") : nothing
S.C[6]>0. ? println("A + B -> A + 2B    with rate  $(S.C[6])") : nothing
S.C[7]>0. ? println("B -> 0     with rate  $(S.C[7])") : nothing
println("\nwith initial condition\nA(0) = $(S.x0[1])\nB(0) = $(S.x0[2])")
return
end

#Usage example:
# S = PathMI.initialize1()
## stores in the mutable structure S of type "System" the reaction network defined in initialize1().
## The parameters of the system can be manually modified, for instance
# S.C[5]=1.5 sets the rate constant of the reaction of index 5 equal to 1.5
# S.x0[1]=20 sets the initial condition of species A at value 20
function initialize1()
    S=create_network()
    S.C[1]=1.   #1# 0->A
    S.C[3]=0.1  #3# A->0
    S.C[5]=1.   #5# A->A+B
    S.C[7]=0.1  #7# B->0
    S.x0=[10;0]
return S
end

function initialize2()
    S=create_network()
    S.C[3]=0.1  #3# A->0
    S.C[5]=1.   #5# A->A+B
    S.C[7]=0.1  #7# B->0
    S.x0=[100;0]
return S
end

function initialize3()
    S=create_network()
    S.C[2]=1.       #1# A->2A
    S.C[4]=0.01     #4# A+B->B
    S.C[6]=0.0002   #6# A+B->A+2B
    S.C[7]=0.1      #7# B->0
    S.x0=[500;100]
return S
end

# OBJECT TO STORE ALL INFORMATION FOR A PATH INSTANCE
mutable struct Path
   status::Bool                 # validity of Path (true=enough allocation size, false=overflow)
   index::Vector{Int64}         # index of current firing channel
   t::Vector{Float64}           # firing times
   x::Matrix{Int64}             # updated state of system after jump
   MA::Matrix{Float64}          # conditional moments of A
   MB::Matrix{Float64}          # conditional moments of B
   pmi::Vector{Float64}         # value of path mutual information sample
   function Path(t0::Float64,x0::Vector{Int64},maxiters::Int64)
        status=true
        index=zeros(Int64,maxiters)
        t=zeros(maxiters)
        t[1]=t0
        x=zeros(Int64,length(x0),maxiters)
        x[:,1]=x0
        MA=zeros(2,maxiters)
        MA[:,1]=[x0[1];x0[1]^2]
        MB=zeros(2,maxiters)
        MB[:,1]=[x0[2];x0[2]^2]
        pmi=zeros(maxiters)
        return new(status,index,t,x,MA,MB,pmi)
   end
end

# Store current state in Path object
function save_path(count::Int64,Xt::Path,index::Int64,time::Float64,xt::Vector{Int64},
                    MomA::Vector{Float64},MomB::Vector{Float64},Infot::Float64)
   Xt.index[count]=index
   Xt.t[count]=time
   Xt.x[:,count]=xt
   Xt.MA[:,count]=MomA
   Xt.MB[:,count]=MomB
   Xt.pmi[count]=Infot
   return
end

# Delete excess allocated elements
function cut_path(Xt::Path, count::Int64)
if count<length(Xt.t)
    Xt.index=Xt.index[1:count]
    Xt.t=Xt.t[1:count]
    Xt.x=Xt.x[:,1:count]
    Xt.MA=Xt.MA[:,1:count]
    Xt.MB=Xt.MB[:,1:count]
    Xt.pmi=Xt.pmi[1:count]
end
return Xt
end
