# This code implements adaptive-bandwidth kernel density estimation and is intended for 
# use on high-throughput sequence data. It accompanies the paper "Adaptive Bandwidth
# Kernel Density Estimation for Next-Generation Sequencing Data" by P. Ramachandran and
# T. J. Perkins (submitted, 2013). This code is free to use for non-commercial purposes.
#
# USAGE:
# 
# In the future, we may organize this code in a more proper package for distribution.
# At present, to use this code, you need merely 'source' this single file. That will
# create five functions: 
#   AKDE which computes an adaptive-bandwidth smooth
#   AssignBandwidths which is called by AKDE to assign a bandwidth to each data point
#   AdaptiveKDE which is called by AKDE to compute the smooth given the bandwidths
#   AKDE_SquareInts computes a square-kernel adaptive KDE smooth, returning the
#      piecewise constant intervals
#   AdaptiveKDE_SquareInts which actually implements the algorithm
#
# If your data is from multiple chromosomes, be aware that this code produces the smooth
# for only one chromosome at a time. Thus, you must call the code with your data from
# each chromosome separately. (We suggest using a for loop in a simple script to loop through the chromosomes!)
#
# The most general function, and possibly the only function you need is AKDE, which
# you use as follows: AKDE(Data,BWRule,BWParam,KernelRule,Lo,Hi)
#   Data is a vector of read positions. They should be positive integers and they should
#      be sorted in increasing order. They are assumed all to be on the same chromosome.
#   BWRule is 0, 1, 2 or 3, and is used to determine the bandwidths associated to each
#      data point. 0 means the same, fixed bandwidth for every point. 1 means the 
#      k-nearest neighbor rule -- the bandwidth associated to each point is equal to 
#      the distance to the kth nearest data point. 2 means the greatest of the 
#      k-nearest neighbor and the distance to the 1-nearest neighbor on either side.
#      (This is different from k-nearest neighbor only if all k-nearest points are on
#      the same side, and the nearest point on the other side is even farther.) 3 means
#      the greater of the kth point to the left and the kth point to the right.
#   BWParam is the bandwidth 'parameter'. For fixed bandwidth (BWRule==0) it is simply
#      the bandwidth. For BWRule 1, 2 or 3, it is the k in the k-nearest neighbor.
#   KernelRule is 1, 2 or 3. 1 is for Gaussian kernels. 2 is for triangle kernels. 
#      3 is for square kernels.
#   Lo and Hi specify the range of basepairs for which the smooth is computed -- it is
#      computed on the range Lo:Hi.
#
# The function AKDE returns the smoothed density over that range. If you want to use
# the smooth for visualization in a genome browser, you may want to output it, say in 
# .wiggle format or the like, which can be done using the Bioconductor package.
# We warn the user that the computations can take some time. A simple fixed-bandwidth
# smooth with square kernels may only take a few minutes to compute on a laptop. However
# using the adaptive KDE rules (which can result in some quite large bandwidths), 
# Gaussian kernels (which also require more computation), or specifying a particularly large
# BWParam, can result in significant computation times. We have benchmarked a
# computation of an adaptive-bandwidth smooth with BWRule==3 BWParam==5 and 
# KernelRule==1 on about 350k reads on Human chromosome 1 (which is ~250Mbp) at about
# one hour on a decently-new Sunfire compute node. You may also need a fair bit of RAM.
# Of course, if you have cluster access, then each chromosome's smooth can be computed 
# on a separate core. This amount of computational effort typically pales in comparison 
# with the effort of generating the data in the first place. Still, if this is too slow
# for you, you could consider our Matlab implementation, which we've benchmarked as 
# running about 10 times faster.
#
# If you are specifically interested in square-kernel smoothing, then another function
# we offer, AKDE_SquareInts, may also be of interest.  For the square kernel, the
# "smooth" is piecewise constant. (We say smooth in quotes because a piecewise constant
# function is not, of course, smooth in the usual continuous or differentiable sense.)
# In any case, this function explicitly computes the intervals of piecewise constant
# height. Usage is: AKDE_SquareInts(Data,BWRule,BWParam,Lo,Hi). The parameters all have
# the same meanings as for the AKDE function. However, the output that is returned is an
# Nx3 matrix where N is the number of piecewise constant components, and the three columns
# give the start position, end position, and height of each component. Zero-height 
# pieces are ommitted from the returned matrix.
#
# We also warn that there is _some_ error checking built into this code, but it is not
# thorough! If you use these functions contrary to how they are intended, you may get
# erroneous/unexpected results without warning.


############
### AKDE ###
############

AKDE <- function(Data,BWRule,BWParam,KernelRule,Lo,Hi) {
  BW <- AssignBandwidths(Data,BWRule,BWParam)
  Smooth <- AdaptiveKDE(Data,BW,KernelRule,Lo,Hi)
  return(Smooth)
}


#######################
### AKDE_SQUAREINTS ###
#######################

AKDE_SquareInts <- function(Data,BWRule,BWParam,Lo,Hi) {
  BW <- AssignBandwidths(Data,BWRule,BWParam)
  Smooth <- AdaptiveKDE_SquareInts(Data,BW,Lo,Hi)
  return(Smooth)
}


########################
### ASSIGNBANDWIDTHS ###
########################

AssignBandwidths <- function(Data, Rule, Param) {
  
  # How many data points?
  N <- length(Data)
  
  # Checking for acceptable parameter values
  if( !(is.numeric(Data)) ) {
    print('AssignBandwidths error: Data is not numeric.')
    return(F)
  }
  if( !((Rule==0) | (Rule==1) | (Rule==2) | (Rule==3)) ) {
    print('AssignBandwidths error: Unknown Rule.')
    return(F)
  }
  if( !(round(Param)==Param) ) {
    print('AssignBandwidths error: Param must be integer.')
    return(F)
  }
  if( Param<=0 ) {
    print('AssignBandwidths error: Param must be positive.')
    return(F)
  }
  if( (Rule>1) & (Param>=N) ) {
    print('AssignBandwidths error: For Rules 2 to 4, Param (k) must be smaller than Data is long.')
    return(F)
  }
  
  
  # Rule 0 means fixed bandwidth of Param for every data point
  if (Rule == 0) {
    BW <- rep(Param,N)
  }
  
  # Rule 1 is classic k-nearest neighbor (Param is k)
  if (Rule == 1) {
    BW <- rep(0,N)
    for (i in 1:N) {
      Lo <- max(1,i-Param)
      Hi <- min(N,i+Param)
      Temp <- sort(abs(Data[Lo:Hi]-Data[i]))
      BW[i] <- Temp[min(Param+1,length(Temp))]
    }
  }
  
  # Rule 2 is the max of Rule 1 and the distance to left and right
  if (Rule == 2) {
    BW <- rep(0,N)
    for (i in 1:N) {
      # Rule 1
      Lo <- max(1,i-Param)
      Hi <- min(N,i+Param)
      Temp <- sort(abs(Data[Lo:Hi]-Data[i]))
      Rule1 <- Temp[min(Param+1,length(Temp))]
      
      # Left neighbor
      Lefti <- max(1,i-1)
      Left <- Data[i]-Data[Lefti]
      
      # Right neighbor
      Righti <- min(i+1,N)
      Right <- Data[Righti]-Data[i]
      
      # And finally...
      BW[i] <- max(Rule1,Left,Right)
    }
  }
  
  # Rule 3 is the max of the kth neighbor to left and kth neighbor to right
  if (Rule == 3) {
    BW <- rep(0,N)
    for (i in 1:N) {
      Lo = max(1,i-Param)
      Hi = min(N,i+Param)
      BW[i] <- max(Data[i]-Data[Lo],Data[Hi]-Data[i])
    }
  }
  
  return(BW)
}


###################
### ADAPTIVEKDE ###
###################

AdaptiveKDE <- function(Data,BW,Kernel,Lo,Hi) {
  
  # How many data points?
  N <- length(Data)
  
  # Initialize smooth
  Smooth <- rep(0,Hi-Lo+1)
  
  # Bandwidth multiplier for support
  if (Kernel==1) BWMult <- 5 # Gaussian
  if (Kernel==2) BWMult <- 2.44949 # Triangle (sqrt(6))
  if (Kernel==3) BWMult <- 1.732051 # Square (sqrt(3))
  
  # Kernel function
  if (Kernel==1) KernelFn <- function(X,Center,Bandwidth) { exp(-(X-Center)^2/(2*Bandwidth^2))/(2.506628*Bandwidth) }
  if (Kernel==2) KernelFn <- function(X,Center,Bandwidth) { (1-abs(X-Center)/(floor(2.44949*Bandwidth)))/floor(2.44949*Bandwidth) }
  if (Kernel==3) KernelFn <- function(X,Center,Bandwidth) { 1/(2*floor(1.732051*Bandwidth)+1) }
  
  # Loop over data points
  for (i in 1:N) {
    
    # Center and bandwidth of kernel
    D <- Data[i]
    B <- BW[i]
    
    # Support of kernel
    SupportLo <- max(Lo,ceiling(D-BWMult*B))
    SupportHi <- min(Hi,floor(D+BWMult*B))
    
    if ( SupportHi >=  SupportLo ) {
      # Explicit support
      Support <- SupportLo:SupportHi
      SupportI <- Support-Lo+1
      
      # Evaluating kernel
      KernelVal <- KernelFn(Support,D,B)
      
      # Adding to smooth
      Smooth[SupportI] <- Smooth[SupportI] + KernelVal
    }
  }
  
  # Normalize
  Smooth <- Smooth/N
  
  return(Smooth)
}


##############################
### ADAPTIVEKDE_SQUAREINTS ###
##############################

AdaptiveKDE_SquareInts <- function(Data,BW,Lo,Hi) {
  # Compute the radius of the square kernel at each data point
  Radii <- floor(sqrt(3)*BW)
  
  # Find start and end positions
  Starts <- pmax(Data-Radii,Lo)
  Ends <- pmin(Data+Radii,Hi)
  
  # We only do ones that overlap the Lo:Hi interval
  # We add an NA to the end of each list, for convenience in the
  # looping logic that follows
  ToDo <- Starts<Ends
  Starts <- Starts[ToDo]
  Ends <- Ends[ToDo]
  Heights <- 1/(2*Radii[ToDo]+1)
  
  # First, make a unified list of starts and ends, and the changes in height
  StartsAndEnds <- c(Starts,Ends)
  HeightChanges <- c(Heights,-Heights)
  
  # Next, sort increasing and add a dummy at the end
  O = order(StartsAndEnds)
  StartsAndEnds <- c(StartsAndEnds[O],Inf)
  HeightChanges <- HeightChanges[O]
  
  # Intervals will be filled in here - start position, end position, height
  Ints <- matrix(NA,nrow=2*length(ToDo),ncol=3)
  NInts <- 0
  CurrIdx <- 1
  CurrHeight <- 0
  while (CurrIdx<length(StartsAndEnds)) {
    # Create the next interval
    IntStart <- StartsAndEnds[CurrIdx]
    
    # Find where interval will end
    NextIdx <- CurrIdx
    while (StartsAndEnds[NextIdx]==StartsAndEnds[CurrIdx]) {
      NextIdx <- NextIdx+1
    }
    
    # As long as we have a real end...
    if ( StartsAndEnds[NextIdx]<Inf ) {
      # Accumlate the new height
      for (i in CurrIdx:(NextIdx-1)) {
        CurrHeight <- CurrHeight + HeightChanges[i]/length(Data)
      }
      
      # Install the next interval, if non-zero
      if (CurrHeight > 0) {
        NInts <- NInts+1
        Ints[NInts,1] <- StartsAndEnds[CurrIdx]
        Ints[NInts,2] <- StartsAndEnds[NextIdx]
        Ints[NInts,3] <- CurrHeight
      }
    }	
    # Update to start of next interval
    CurrIdx <- NextIdx
  }
  
  Ints <- Ints[1:NInts,]
  return(Ints)
}