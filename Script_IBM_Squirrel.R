# Individual-based model for ground squirrel 
# example from section 2 in Coste et al. 2021
# The Kinship Matrix: Inferring the Kinship Structure of a Population From Its Demography

# Sarah Cubaynes - Juan Pablo Ramirez Loza
# June 2021

# Packages used to run the simulations
library(NetLogoR)
library(testthat)
library(kinship2)
library(SciViews)
library(rlist)

# 1- Define the Leslie matrix
# pre breedng: a1 yearling, a2 2yo, a3 adults
nstages = 3
F1 = 0.23
F2 = 0.51
F3 = 0.84
s1 = 0.573
s2 = 0.537
s3 = 0.460
Mat <- matrix(c(F1,F2,F3,
                s1,0,0,
                0,s2,s3), byrow=T, nrow=3)

# 2 -Calculate growth rate and stable age structure with popbio
library(popbio)
lambda_equi <- lambda(Mat)  # asymtptotic growth rate
w_equi <- stable.stage(Mat) # stable age distribution (% of each age class in the pop at equilibrium)

# 3 -Load expected kinship matrices from Coste et al. model
# expected kinship matrices 
expectedK <- list(mother = matrix(c(0,	0,	0,
                                    0.1311,	0,	0,
                                    0.3774,	0.2444,	0.0770),
                                   byrow = T,
                                  nrow=3),
                  gmother =matrix(c(0,	0,	0,
                                    0,	0,	0,
                                    0.1059,    0.0487,    0.0153),
                                  byrow = T,
                                  nrow=3),
                  sister =matrix(c(  0.6037,	0.3839,	0.1406,
                                     0.2189,	0.3441,	0.1550,
                                     0.0802,	0.1551,	0.2250),
                                 byrow = T,
                                 nrow=3),
                  aunts = matrix(c(0.1790,    0.0890,    0.0280,
                                   0.1812,    0.1020,    0.0347,
                                   0.2002,    0.1899,    0.1085),
                                 byrow = T,
                                 nrow=3),
                  children = matrix(c(0,	0.23,	0.6618,
                                      0,	0,	0.243,
                                      0,	0,	0.0770),
                                    byrow = T,
                                    nrow=3),
                  gchildren = matrix(c(0,	0,	0.18586,
                                       0,	0,	0.048972,
                                       0,	0,	0.015337),
                                     byrow = T,
                                     nrow=3), 
                  cousins = t(matrix(c(0.3316,    0.3017,    0.1783,
                                       0.1720,    0.1890,    0.1396,
                                       0.1017,    0.1397,    0.1579),
                                     byrow = T,
                                     nrow=3)),
                  nieces =matrix(c(0.17895,	0.25743,	0.2883,
                                   0.050909,	0.10254,	0.15564,
                                   0.015944,	0.034673,	0.097394),
                                 byrow = T,
                                 nrow=3) )

# 4 - Define functions needed in the IBM (derived from Bauduin et al. 2020)
# Sample function that works the same with 1 (length(x) == 1) or several items (length(x) >= 1) to sample
sample.vec <- function(x, ...) x[sample(length(x), ...)]

# Reproduction
repro <- function(squirrels, allsquirrelsID, yearSim, popSim){
  a1Reproduce <- NLwith(agents = squirrels, var = "stage", val = 1)
  a2Reproduce <- NLwith(agents = squirrels, var = "stage", val = 2)
  a3Reproduce <- NLwith(agents = squirrels, var = "stage", val = 3)
  if(sum(NLcount(a1Reproduce),NLcount(a2Reproduce),NLcount(a3Reproduce)) != 0){
    if(runTests){
      numsquirrels <- NLcount(squirrels)
    }
    IDa1Reproduce <- of(agents = a1Reproduce, var = "who")
    IDa2Reproduce <- of(agents = a2Reproduce, var = "who")
    IDa3Reproduce <- of(agents = a3Reproduce, var = "who")
    # Number of pups each female will have (different number per female)
    nPups1 <-  rpois(n = length(IDa1Reproduce), lambda = F1)
    nPups2 <-  rpois(n = length(IDa2Reproduce), lambda = F2)
    nPups3 <-  rpois(n = length(IDa3Reproduce), lambda = F3)
    nPups <- c(nPups1 , nPups2 , nPups3)
    if(sum(nPups) != 0){
      # Remove from the loop the females which have 0 pup
      IDa1Reproduce <- IDa1Reproduce[nPups1 != 0] 
      nPups1 <- nPups1[nPups1 != 0]
      a1Reproduce <- turtle(turtles = squirrels, who = IDa1Reproduce)

      IDa2Reproduce <- IDa2Reproduce[nPups2 != 0] 
      nPups2 <- nPups2[nPups2 != 0]
      a2Reproduce <- turtle(turtles = squirrels, who = IDa2Reproduce)
      
      IDa3Reproduce <- IDa3Reproduce[nPups3 != 0] 
      nPups3 <- nPups3[nPups3 != 0]
      a3Reproduce <- turtle(turtles = squirrels, who = IDa3Reproduce)
      
      IDFemaleReproduce <- c(IDa1Reproduce,IDa2Reproduce,IDa3Reproduce)
      nPups <- c(nPups1, nPups2, nPups3)
      
      # Create the new squirrels with hatch()
      squirrels <- hatch(turtles = squirrels, who = IDFemaleReproduce, n = nPups, breed = "newborn") # breed = "newborn" to recognize the pups in the squirrels object
      newborn <- NLwith(agents = squirrels, var = "breed", val = "newborn")
      # The newborns inherit all the parent's (femaleReproduce) parameters except for the who numbers
      # Some of the inherited variables must be changed
      # The who numbers (IDs) also need to be updated so that newborns never have an IDs of an already dead wolf from this population (this is needed to later define the pedigree)
      uniqueWho <- seq(from = max(allsquirrelsID) + 1, to = max(allsquirrelsID) + NLcount(newborn), by = 1)
      
      # Update the new squirrels variables
      squirrels <- NLset(turtles = squirrels, agents = newborn, 
                      var = c("stage", "breed", "who", "motherID", "cohort"),
                      val = data.frame(stage = 0,
                                       breed = "turtle",
                                       who = uniqueWho,
                                       motherID = rep(IDFemaleReproduce, nPups),
                                       cohort = yearSim))
      
      # Update allsquirrelsID with the new squirrels ID
      allsquirrelsID <- c(allsquirrelsID, uniqueWho)
      
      if(runTests){
        # Make sure pups were integrated in the squirrels object
        expect_equivalent(NLcount(squirrels), numsquirrels + sum(nPups)) 
      }
    }
  }
  
  if(runTests){
    # There should not be duplicated IDs
    expect_equal(length(allsquirrelsID), length(unique(allsquirrelsID)))
  }
  
  temporaryOutputs <- popSim
  temporaryOutputs[[length(popSim) + 1]] <- squirrels
  
  return(list(squirrels, allsquirrelsID))
}


# Mortality
mortality <- function(squirrels){
  if(runTests){
    numsquirrels <- NLcount(squirrels) 
  }
  # Pup Mortality is fixed to 0 here, because pup mortality is included in the fecundity parameter
  # so all pups produced at reproduction year t become the yearling a1 at the beginning of next time interval
  
  # Yearling mortality
  squirrelsYearling <- NLwith(agents = squirrels, var = "stage", val = 1) # because aging process occurs before mortality, yearling mortality applies to individuals aged 2
  IDYearling <- of(agents = squirrelsYearling, var = "who")
  probDeathYearling <- 1 -s1 #rnorm(1, mean = mortalityYearling, sd = mortalityYearlingSD)
  if(probDeathYearling < 0){probDeathYearling <- mortalityYearling} # to avoid negative probability
  deadYearling <- rbinom(n = length(IDYearling), size = 1,
                         prob = probDeathYearling) # 0 = survive, 1 = die
  squirrels <- die(turtles = squirrels, who = IDYearling[deadYearling == 1])
  
  # Two year old mortality
  squirrels2Y <- NLwith(agents = squirrels, var = "stage", val = 2) # because aging process occurs before mortality, yearling mortality applies to individuals aged 2
  ID2Y <- of(agents = squirrels2Y, var = "who")
  probDeath2Y <- 1 -s2 #rnorm(1, mean = mortalityYearling, sd = mortalityYearlingSD)
  if(probDeath2Y < 0){probDeath2Y <- mortality2Y} # to avoid negative probability
  dead2Y <- rbinom(n = length(ID2Y), size = 1,
                         prob = probDeath2Y) # 0 = survive, 1 = die
  squirrels <- die(turtles = squirrels, who = ID2Y[dead2Y == 1])

  # Older adult mortality
  squirrels3Y <- NLwith(agents = squirrels, var = "stage", val = 3) # because aging process occurs before mortality, yearling mortality applies to individuals aged 2
  ID3Y <- of(agents = squirrels3Y, var = "who")
  probDeath3Y <- 1 -s3 #rnorm(1, mean = mortalityYearling, sd = mortalityYearlingSD)
  if(probDeath3Y < 0){probDeath3Y <- mortality3Y} # to avoid negative probability
  dead3Y <- rbinom(n = length(ID3Y), size = 1,
                   prob = probDeath3Y) # 0 = survive, 1 = die
  squirrels <- die(turtles = squirrels, who = ID3Y[dead3Y == 1])
  
  if(runTests){
    # All dead individuals were removed from squirrels
    expect_equal(NLcount(squirrels),
                 numsquirrels - NLcount(deadYearling) - sum(dead2Y) - sum(dead3Y) )
  }
  return(squirrels)
}

# Aging
aging <- function(squirrels){
  # stage transition (squirrels in stage 3 remain in stage 3)
  transients <- NLwith(agents=squirrels,
               var="stage",
               val=0:2)
  stagesquirrel <- of(agents = transients, var = "stage")

  squirrels <- NLset(turtles = squirrels, 
                     agents=transients , 
                     var= "stage", 
                     val = stagesquirrel+1)
  if(runTests){
    # No squirrel should be in stage > 3
    expect_true(all(stagesquirrel <= nstages))
  }
  
  return(squirrels)
}

# Calculate Kinship structure at the last time step of the simulated population 
KinshipStructure <- function(popSim, nstages){ # entry = popSim which contains the result of one simulation replicate from the IBM
  # popSim contains all indiv that ever lived in the pop, it is a list with t elements (pop at each time step)
  # For individuals ever alive in the pop (dead or not at last time step)
  dataall  <- unique(do.call(rbind, popSim)  ) # all indiv ever lived in population
  allwithparents <- dataall[!is.na(dataall$motherID),] #se data for individuals with known parents only
  # For individuals alive at last time step
    datalast <- as.matrix(popSim[[length(popSim)]])
    if(nrow(datalast)==0){
      kinshipinfo <- NA
    }else{
  if(is.null(popSim$Sex)==TRUE){
    AllFem <- allwithparents
    datalastF <- datalast[ (!is.na(datalast$motherID) ) ,] # keep only females with known mother 
    
  }else{
    # keep only female with know mother ID (all except founder individuals)
    AllFem <- allwithparents[allwithparents$sex=="F",] 
    datalastF <- datalast[ (datalast$sex=="F" & !is.na(datalast$motherID) ) ,] # keep only females with known mother 
  }
  # Ages - calculated among all females with known mother ever alive
  stages<- 1:nstages
  # ID of all mothers and grandmothers ever lived in pop
  IDMothers <- unique(AllFem[,'motherID'] ) # mothers of at least one female
  IDGrandMothers <- unique( AllFem[which(AllFem$who %in% IDMothers),'motherID'] )# grand mother of at least one female
  # Mothers alive at last time step
  IDMothersAlive <- datalastF$who[datalastF$who %in% IDMothers]
  # Ages of mothers alive at last time step
  
  # Grand Mothers alive at last time step
  IDGrandMothersAlive <- datalastF$who[datalastF$who %in% IDGrandMothers]
  # Create objects to store quantitities  
  #Mother & Grandmother
  nb_MotherAlive <- kinship_Mother <- matrix(NA, nrow=max(stages),ncol=max(stages))
  nb_GrandmotherAlive <- kinship_Grandmother <- matrix(NA, nrow=max(stages),ncol=max(stages))
  n_MA <- n_GMA <-Pr_MA_egoagei <- Pr_GMA_egoagei <-NULL
  mums <- list()
  #Sisters
  kinship_Sisters <- matrix(NA, nrow=max(stages),ncol=max(stages))
  #Aunts
  kinship_Aunts <- matrix(NA, nrow=max(stages),ncol=max(stages))
  # Children
  kinship_Children <- matrix(NA, nrow=max(stages),ncol=max(stages))
  # Grand Children
  kinship_GrandChildren <- matrix(NA, nrow=max(stages),ncol=max(stages))
  # Cousins
  kinship_Cousins <- matrix(NA, nrow=max(stages),ncol=max(stages))
  #Nieces
  kinship_Nieces <- matrix(NA, nrow=max(stages),ncol=max(stages))
  #Summary of results
  Res_summary <- data.frame(ego.age=stages,
                            nb.ego=c(rep(NA,length(stages))),  # number of ego aged i at last time step
                            nb.ego.MA=c(rep(NA,length(stages))),  # number of ego aged i with mother alive at last time step
                            Pr_MA =c(rep(NA,length(stages))), # proba mother alive for ego age i
                            nb.ego.GMA=c(rep(NA,length(stages))),  # number of ego aged i with mother alive at last time step
                            Pr_GMA =c(rep(NA,length(stages))),
                            Avg_Sis = c(rep(NA,length(stages))),
                            Avg_Aunts = c(rep(NA,length(stages))),
                            Avg_Ch = c(rep(NA,length(stages))),
                            Avg_GCh = c(rep(NA,length(stages))),
                            Avg_Cou = c(rep(NA,length(stages))),
                            Avg_Nie = c(rep(NA,length(stages)))
  )
  
  # Calculate kinship structure for ego age i, each kin age j
  for(i in stages){  #loop over ego age
    ego_agei <- datalastF[datalastF$stage==i,] # ego age i alive at last time step
    
    if(nrow(ego_agei)==0){
      nb_MotherAlive[i,] <- rep(NA,length(stages))
      nb_GrandmotherAlive[i,] <-rep(NA,length(stages))
      kinship_Mother[i,] <- rep(NA, length(stages))
      kinship_Grandmother[i,] <- rep(NA, length(stages))
      kinship_Sisters[i,] <- rep(NA,length(stages))
      kinship_Aunts[i,] <- rep(NA,length(stages))
      kinship_Children[i,] <- rep(NA,length(stages))
      kinship_GrandChildren[i,] <- rep(NA,length(stages))
      kinship_Cousins[i,] <- rep(NA,length(stages))
      kinship_Nieces[i,] <- rep(NA, length(stages))
      
      Res_summary[i,"nb.ego"] <- nrow(ego_agei)
      Res_summary[i,"nb.ego.MA"] <- NA
      Res_summary[i,"Pr_MA"] <- NA
      Res_summary[i,"nb.ego.GMA"] <- NA
      Res_summary[i,"Pr_GMA"] <-NA
      Res_summary[i,"Avg_Sis"] <- NA
      Res_summary[i,"Avg_Aunts"] <- NA
      Res_summary[i,"Avg_Ch"] <- NA
      Res_summary[i,"Avg_GCh"] <- NA
      Res_summary[i,"Avg_Cou"] <- NA
      Res_summary[i,"Avg_Nie"] <- NA
      
    }else if(nrow(ego_agei)>0){ 
      Mothers.ID.ego <- ego_agei$motherID # mothers ID of ego age i
      Mothers.Alive.ID <- ego_agei$motherID[which(Mothers.ID.ego %in% IDMothersAlive)] # alive mothers ID of ego age i, repeated by number of daughters
      mums[[i]] <- Mothers.Alive.ID
      n_MA[i] <- length(Mothers.Alive.ID) # nb of ego age i with mother alive at last time step
      Pr_MA_egoagei[i]  <- n_MA[i]/ nrow(ego_agei)  # proba mother alive for ego age i
      
      GM.ID <- AllFem$motherID[match(Mothers.ID.ego, AllFem$who)] # GM ID of ego age i
      GM.Alive.ID <- GM.ID[GM.ID %in% IDGrandMothersAlive]
      n_GMA[i] <- length(GM.Alive.ID) # nb of ego age i with grandmother alive at last time step
      Pr_GMA_egoagei[i]  <- n_GMA[i]/ nrow(ego_agei)  # proba mother alive for ego age i
      
      # age of mother AND GRAND MOTHER if alive at last time step
      Mage <- datalastF$stage[match(Mothers.Alive.ID, datalastF$who )] # match ID of mothers alive at last time step in datalast (which includes all females alive at last time step)
      GMage <- datalastF$stage[match(GM.Alive.ID, datalastF$who )] # match ID of mothers alive at last time step in datalast (which includes all females alive at last time step)
      
      for(j in stages){ # loop over ages
        ind_agej <- datalastF[datalastF$stage==j] #individuals of age j alive at last time step
        
        # FOR MOTHERS
        nb_MotherAlive[i,j] <- sum(Mage==j) 
        kinship_Mother[i,j] <- sum(Mage==j) / nrow(ego_agei) # kinship matrix Pr of mother alive
        # FOR GRAND MOTHERS
        nb_GrandmotherAlive[i,j] <- sum(GMage==j) 
        kinship_Grandmother[i,j] <- sum(GMage==j) / nrow(ego_agei) # kinship matrix Pr of mother alive
        
        
        GM.IDs.ind_agej <-AllFem$motherID[match(ind_agej$motherID, AllFem$who )] # GM ID of ind age j alive at last time step, the GM can can be dead or alive
        nb.aunts.agej.ego <- nb.sis.agej.ego <- nb.children.agej.ego <- nb.grandchildren.agej.ego <- nb.cousins.agej.ego <- nb.nieces.agej.ego <- NULL # empty object to store results
        
        for(e in 1:nrow(ego_agei)){ # loop over each ego aged i
          
          ego.GM.ID <- GM.ID[e] #ego's grand mother ID, dead or alive
          
          #For nieces
          DaugthersOfEgoMother <- AllFem[which(AllFem$motherID==ego_agei$motherID[e]),] #Get all the daughters of ego's mother
          AllSistersDeadorAlive <- DaugthersOfEgoMother[-which(DaugthersOfEgoMother$who==ego_agei$who[e]),] #get only ego's sisters (eliminate ego)
          nb.nieces.agej.ego[e] <- sum(ind_agej$motherID %in% AllSistersDeadorAlive$who)
          
          
          nb.children.agej.ego[e] <- sum(ind_agej$motherID %in% ego_agei$who[e]) # nb of children of ego e that are age j
          nb.grandchildren.agej.ego[e] <- sum(GM.IDs.ind_agej %in% ego_agei$who[e]) # nb of grandchildren of ego e that are age j
          
          if (sum(is.na(ego.GM.ID))==1) { #Added
            nb.aunts.agej.ego[e] <- 0       #Added
          } else if (sum(is.na(ego.GM.ID))==0) {   #Added
            
            if (sum(ind_agej$who==Mothers.ID.ego[e])==0){ # if ego's mother is dead
              nb.aunts.agej.ego[e] <- sum(ind_agej$motherID == GM.ID[e] )
            } else if (sum(ind_agej$who==Mothers.ID.ego[e])>0){ # if ego's mother is alive
              ind_agej_withoutegosmother <- ind_agej[-which(ind_agej$who==Mothers.ID.ego[e]),] # all indiv alive at last time step except for ego's mother
              nb.aunts.agej.ego[e] <- sum(ind_agej_withoutegosmother$motherID == GM.ID[e])
            }  #Added
          }
          
          if (j==i){ #same cohort sisters
            nb.sis.agej.ego[e] <- sum(ind_agej$motherID==Mothers.ID.ego[e]) - 1
            
            if (sum(is.na(ego.GM.ID))==1) { #Added
              nb.cousins.agej.ego[e] <- 0         #Added
            } else if (sum(is.na(ego.GM.ID))==0) {      #Added
              
              nb.cousins.agej.ego[e] <-sum(GM.IDs.ind_agej %in% ego.GM.ID) - nb.sis.agej.ego[e] - 1 # nb of cousins of ego e that are age j
            }  #Added
          } else if (j!=i){ 
            nb.sis.agej.ego[e] <- sum(ind_agej$motherID==Mothers.ID.ego[e])
            
            if (sum(is.na(ego.GM.ID))==1) {     #Added
              nb.cousins.agej.ego[e] <- 0        #Added
            } else if (sum(is.na(ego.GM.ID))==0) {     #Added
              
              nb.cousins.agej.ego[e] <-sum(GM.IDs.ind_agej %in% ego.GM.ID) - nb.sis.agej.ego[e]
            }
          }  #Added
        } # end loop on ego e
        
        kinship_Sisters[i,j] <- mean(nb.sis.agej.ego)
        kinship_Aunts[i,j] <- mean(nb.aunts.agej.ego)
        kinship_Children[i,j] <- mean(nb.children.agej.ego)
        kinship_GrandChildren[i,j] <- mean(nb.grandchildren.agej.ego)   
        kinship_Cousins[i,j] <-  mean(nb.cousins.agej.ego)
        kinship_Nieces[i,j] <- mean(nb.nieces.agej.ego)
      } # end loop on sister's and aunt's ages
      Res_summary[i,"nb.ego"] <- nrow(ego_agei)
      Res_summary[i,"nb.ego.MA"] <- n_MA[i]
      Res_summary[i,"Pr_MA"] <- Pr_MA_egoagei[i]
      Res_summary[i,"nb.ego.GMA"] <- n_GMA[i]
      Res_summary[i,"Pr_GMA"] <- Pr_GMA_egoagei[i]
      Res_summary[i,"Avg_Sis"] <- sum(kinship_Sisters[i,])
      Res_summary[i,"Avg_Aunts"] <- sum(kinship_Aunts[i,])
      Res_summary[i,"Avg_Ch"] <- sum(kinship_Children[i,])
      Res_summary[i,"Avg_GCh"] <- sum(kinship_GrandChildren[i,])
      Res_summary[i,"Avg_Cou"] <- sum(kinship_Cousins[i,])
      Res_summary[i,"Avg_Nie"] <- sum(kinship_Nieces[i,])
    } # end else if
    
  } # end loop on ego's age
  
  rsummary <- round(Res_summary,2) # summary per ego age
  
  kinshipinfo <- list(res=rsummary,
                      K_Mothers=kinship_Mother,
                      K_GrandMothers=kinship_Grandmother,
                      K_Sisters=kinship_Sisters,
                      K_Aunts=kinship_Aunts,
                      K_Children=kinship_Children,
                      K_GrandChildren=kinship_GrandChildren,
                      K_Cousins=kinship_Cousins, #K_Coousins
                      K_Nieces=kinship_Nieces)
  }
  return(kinshipinfo)
}



# 5- Run the IBM simulations

# Define the model parameters
nReplicate <- 300 # how many replicates of the population simulation
nYearSim <- 50 # how many years of simulation

# Create the initial population at the beginning of each simulation replicate
# Create a data frame of the individual wolf characteristics in the initial population
initPop_DF <- data.frame(ID = 1:100, # 100 individuals with id from 1 to 100
                             stage = c(rep(c(3, 2, 1, 3, 2), 20))
)

# Create the initial population using the NetLogoR package
# Create an agentMatrix object
# Individual IDs are "who" in the agentMatrix and starts at 0 (automatically created when creating the individuals)
init <- function(initPop_DF){
  
  # Initialize the model objects (i.e., land and squirrels)
  land <- createWorld(0, 100, 0, 100) # create a fictive land (does not matter, will not be used)
  squirrels <- createTurtles(n = nrow(initPop_DF), world = land) # create as many squirrels as nrow(initPopWolf_DF)
  squirrels <- turtlesOwn(turtles = squirrels, tVar = "stage", tVal = initPop_DF[, "stage"])
  squirrels <- turtlesOwn(turtles = squirrels, tVar = "motherID", tVal = as.numeric(NA))
  squirrels <- turtlesOwn(turtles = squirrels, tVar = "cohort", tVal = as.numeric(NA))
  return(squirrels)
}

# Running tests during the simulation to identify potential bugs
runTests <- FALSE # put FALSE for faster simulations

# Create the output files
popSimRep <- kinship <- lambda <- lambda_lasttimestep <- w <- list() # record as the list the state of the population each year for each simulation
for(j in 1:nReplicate){ # loop over the simulation replicates
  
  squirrels <- init(initPop_DF)
  
  # Create a list to record the state of the population each year (for THIS replicate)
  popSim <- list()
  popSim[[1]] <- squirrels # record the initial state of the wolf population
  
  # Create some vectors needed for the sub-models
  allsquirrelsID <- of(agents = squirrels, var = "who") # keep in memory all the squirrels ID ever created
  yearSim <- 0 # used to indicate the cohorts, the number itself doesn't matter, it just needs to be updated at each loop
  
  for(i in 1:nYearSim){ # loop over the years simulated
    
    # Update the year simulated at the beginning of each year
    yearSim <- yearSim + 1
    if(NLcount(squirrels) != 0){ # do not execute the sub-models if there are no squirrels alive
      
      ##############
      # Reproduction
      resRepro <- repro(squirrels, allsquirrelsID, yearSim, popSim) # result = list(squirrels, allsquirrelsID, allsquirrelsRelatedness). yearSim  and popSim are used but not modified so not returned
      squirrels <- resRepro[[1]] # update the objects each time there are modified
      allsquirrelsID <- resRepro[[2]]
      ##############
      
       ###########
      # Mortality
      squirrels <- mortality(squirrels)
      ###########
 
      #######
      # Aging
      squirrels <- aging(squirrels)
      #######
      
      # Add the new state of the population
      popSim[[i + 1]] <- squirrels
      
    } else { # if there are no squirrels anymore
      # Do not run the sub-models but still add the new state of the population (i.e., no squirrels)
      popSim[[i + 1]] <- squirrels
    }
    
  } # end of the year simulated
  
  
  # Add the list popSim to the list popSimRep
  popSimRep[[j]] <- popSim # population simulated over all the years for the current simulation replicate
  # Print the number of the current replicate
  print(paste0("Replicate", j))
  
  # At the end of each replicate, save popSimRep and packDyn in a .Rdata file so that if there a problem, you can access to all simulation replicates already done
  #save(popSimRep, file = "wolfSimRes.RData")
  
  # Calculate kinship matrices
  kinship[[j]] <- KinshipStructure(popSim, nstages)
  
  # At the end of each replicate, save popSimRep and packDyn in a .Rdata file so that if there a problem, you can access to all simulation replicates already done
  #save(kinship, file = "kinship_wolfSimRes.RData")
  
  # growth rate
  N <-unlist(lapply(popSim,nrow))
  lambda[[j]] <- N[2:nYearSim] / N[1:(nYearSim-1)]
  lambda_lasttimestep[j] <- N[nYearSim] / N[(nYearSim-1)]
  
  # stable age structure
  if(nrow(popSim[[nYearSim]])!=0){
  w[[j]] <- as.data.frame(table(popSim[[nYearSim]]$stage))[,2] /N[nYearSim]
  }else{
    w[[j]] <- c(0,0,0)  
  }
} # end of the replicate simulation replicate
###########################

# SAVE RESULTS
save(kinship, file = "kinship_squirrel_T50.RData")
save(popSim,file="popSim_squirrel_T50.Rdata")
save(lambda,file="lambda_squirrel_T50.Rdata")
save(w,file="StableAge_squirrel_T50.Rdata")

# Check pop sizes - select only replicates with enough sample to calculate kinship (>100 indiv at the last time step)
poplast <- lapply(popSimRep, tail, n = 1L)
poplastmat <- sapply(poplast,as.matrix)
nlastmat <- sapply(poplastmat,nrow)
hist(nlastmat)
more100ind <- which(nlastmat>50)[1:200]

## ANALYSE RESULTS

# check stable age structure
w1 <-list.rbind(w)[,1]# quantile(,probs=c(0.025,0.5,0.975),na.rm=T)
w2 <-list.rbind(w)[,2]#quantile(,probs=c(0.025,0.5,0.975),na.rm=T)
w3 <-list.rbind(w)[,3]#quantile(,probs=c(0.025,0.5,0.975),na.rm=T)
summary(w1)
summary(w2)
summary(w3)
w_equi
boxplot(w1,w2,w3,names=c("a1","a2","a3"),
        xlab="Stage",
        ylab="Proportion of individuals",
        las=TRUE) # simulated value for each age class
points(1:nstages,w_equi,col="red",pch=15,cex=1) # expected value in red

# check lambda at each time step to make sure that the number of time of steps is enough to reach convergence
lamda_lastT <-list.rbind(lambda)[,(nYearSim-1)]
summary(lamda_lastT)
lambda_equi  # expected lambda
# View lambda calculated for each replicate and expected lambda in red
plot(density(lamda_lastT,na.rm=T),
     main="",
     xlab="Growth rate",
     las=T)
abline(v=lambda_equi,lwd=2,col="red")


# Select kinship for 200 replicates for which population did not go extinct
notextinct <- which(sapply(w,sum)!=0)[1:200]
nReplicate <- length(notextinct)
kinship<- kinship[notextinct]

# Get mean/median with 95%CI simulated kinship matrices 
nkins <- 8 # mother, grandmother, sister, aunts, children, grandchildren, 
#cousins, nieces
summarymat <- Kmat <- simulatedK <- quantK <- list()

for(j in 1:nkins){
    for (i in 1:nReplicate){
  summarymat[[i]] <- t(as.matrix(kinship[[i]][[1]]) )
    Kmat[[i]] <- t(as.matrix(kinship[[i]][[j+1]]) )
    }
  simulatedK[[j]] <- round(mean.list(Kmat,na.rm=TRUE),3) # average kinship matrix over replicates
  quantK[[j]] <- apply(simplify2array(Kmat), 1:2, quantile, prob = c(0.025, 0.5,0.975),na.rm=TRUE) # get confidence intervals
 summaryres <- round(mean.list(summarymat,na.rm=TRUE),3)
}
names(simulatedK) <-c("mother", "grandmother", "sisters", "aunts", "children", 
                      "grandchildren","cousins", "nieces") 

# summary of results
# in row : age of ego
# in column : age of kin
simulatedK # list of mean kinship matrix for each kin in this order
#quantK # list of 2.5% , median and 97% quantile for each kin in this order
# mother, gmother, sisters, aunts, female children, female gchildren, female cousins, nieces

# Compare simulated versus expected kinship matrices
mean(round(simulatedK$mother - expectedK$mother ,2))
mean(round(simulatedK$aunts - expectedK$aunts,2))

# Plot simulated versus expected kinship structure for each stage of ego (averaged over kin age)
  Mother_simulated_expected <- data.frame(stage=c("a1","a2","a3"),
                         meanvalue=apply(simulatedK$mother,2,sum),
                         cilow= apply(quantK[[1]][1,,],2,sum),
                         cihigh=apply(quantK[[1]][3,,],2,sum) ,
                         predicted= apply(expectedK$mother,2,sum) )
  p <- ggplot(Mother_simulated_expected, aes(x=stage, y=meanvalue)) + 
    geom_point(size=2, position = position_nudge(x = -0.05))+
    geom_errorbar(aes(ymin=cilow, ymax=cihigh), width=.1, position = position_nudge(x = -0.05))+  
    geom_point(aes(x=stage, y=predicted), size=2, shape=17,col="red", position = position_nudge(x = +0.05)) 
  p +  ylim(0,1) +
    ylab("Probability of mother being alive")
  ggsave("Mother.pdf")

  Aunts_simulated_expected <- data.frame(stage=c("a1","a2","a3"),
                                       meanvalue=apply(simulatedK$aunts,2,sum),
                                       cilow= apply(quantK[[4]][1,,],2,sum),
                                       cihigh=apply(quantK[[4]][3,,],2,sum) ,
                                       predicted= apply(expectedK$aunts,2,sum) )

    p2 <- ggplot(Aunts_simulated_expected, aes(x=stage, y=meanvalue)) + 
    geom_point(size=2, position = position_nudge(x = -0.05))+
    geom_errorbar(aes(ymin=cilow, ymax=cihigh), width=.1, position = position_nudge(x = -0.05))+  
    geom_point(aes(x=stage, y=predicted), size=2, shape=17,col="red", position = position_nudge(x = +0.05)) 
  p2 +  ylim(0,1.5) +
    ylab("Number of aunts alive")

    ggsave("Aunts.pdf")
