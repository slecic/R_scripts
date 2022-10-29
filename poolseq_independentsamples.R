# calculate relatedness for all 2014 samples (814 + 1200)
# read in the files
# 814 pedigreed individuals
samplefulluniqueall <- read.csv("/Users/sonjalecic/Documents/samplefulluniqueall.csv")
# additional 1200 individuals
addindPed <- 
  read.csv("/Users/sonjalecic/Desktop/rereadyfortakeoff/microsatsNonsampleBreedersPed.csv") 

# compile a pedigree file with sample ID, father ID, mother ID and sex
AllSamp <- data.frame(ID = c(samplefulluniqueall$ID, addindPed$ID),
                      father = c(samplefulluniqueall$genfather, addindPed$genfather),
                      mother = c(samplefulluniqueall$genmother, addindPed$genmother),
                      sex = c(samplefulluniqueall$sex, addindPed$sex))

suppressMessages(library(kinship2))
# calculate kinship
AllSamplesKinship <- kinship(AllSamp$ID, AllSamp$father, 
                             AllSamp$mother, AllSamp$sex)


# multiply everything with 2 to get the relatedness values because relatednes = kinship * 2
AllSamplesKinship <- 2*AllSamplesKinship
hist(AllSamplesKinship[AllSamplesKinship > 0 & AllSamplesKinship < 1], 
     breaks = 200,
     xlab = "relatedness",
     ylab = "pairs of individuals",
     ylim = c(0, 25000),
     main = "relatedness of all individuals",
     col = "plum3")

mean(AllSamplesKinship)
max(AllSamplesKinship) # > 1 -> correction for inbreeding is needed
# save the files
write.table(AllSamplesKinship, file="AllSamplesKinship.Rdata", 
            row.names = T, col.names = T)
write.table(AllSamplesKinship, "AllSamplesKinship.txt", 
            row.names = T, col.names = T)

# calculate individual inbreeding coefficient
suppressMessages(library(pedigree))
ped <- subset(AllSamp, select = -sex) # exlude the sex column
ord <- orderPed(ped) # order the ped file
pedigree <- ped[order(ord),]
inbreeding <- calcInbreeding(pedigree)  # calculate inbreeding
indivinbreed <- data.frame(ID=pedigree$ID, F=inbreeding)
table(indivinbreed$F)
mean(indivinbreed$F)

AllSamplesKinshiplist <- 
  melt(AllSamplesKinship) # melt the relatedness matrix into an adjecency list 
AllSamplesKinshiplist$X1F <- NA
AllSamplesKinshiplist$X2F <- NA 
AllSamplesKinshiplist$X1F <- indivinbreed$F[match(AllSamplesKinshiplist$X1, indivinbreed$ID)]
AllSamplesKinshiplist$X2F <- indivinbreed$F[match(AllSamplesKinshiplist$X2, indivinbreed$ID)]
# check
#AllSamplesKinshiplist[AllSamplesKinshiplist$X1F > 0, ]

# correct R accoring to the formula R(x,y)= 2f(x,y)/sqrt((1+Fx)*(1+Fy)) 
# from Galla et al., 2020, Evol. Appl.
AllSamplesKinshiplist$r <- 
  AllSamplesKinshiplist$value/(sqrt((1 + AllSamplesKinshiplist$X1F)*
                                      (1 + AllSamplesKinshiplist$X2F)))

# check the maximum value
max(AllSamplesKinshiplist$r) # it is 1 as it should be (individuals againts themselves)

# save the adjecency list with corrected r values
AllSamplesKinshipCorrected <- AllSamplesKinshiplist[, c(1, 2, 6)] 

# convert the adjecency list into an adjecency matrix
suppressMessages(library(circlize))
allsampAdjmatrix <- 
  adjacencyList2Matrix(AllSamplesKinshipCorrected, square = TRUE) # takes a couple of min
suppressMessages(library(gdata))
lowerTriangle(allsampAdjmatrix) <- upperTriangle(allsampAdjmatrix, byrow=TRUE)
diag(allsampAdjmatrix) <- 1
# check the max and the min relatedness
max(allsampAdjmatrix) # 1
min(allsampAdjmatrix) # 0
write.table(allsampAdjmatrix, file="allsampAdjmatrix.Rdata", row.names = T, col.names = T)
write.table(allsampAdjmatrix, file="allsampAdjmatrix.txt", row.names = T, col.names = T)

# read in 1130 independent individuals (output of fastindep; threshold 0.5)
AllSampIndep <- 
  scan("/Users/sonjalecic/Documents/AllSampIndep.txt",# 1130 = maximum indep individulas
       what = "character") 
length(AllSampIndep)

# extract those 1130 individuals from the original relatedness matrix 
allsampAdjmatrix <-read.table("allsampAdjmatrix.Rdata") 
AllSampIndepMatrix <- allsampAdjmatrix[AllSampIndep, AllSampIndep]
# and check relatedness values (if all values are below 0.5)
dim(AllSampIndepMatrix) # 1130 x 1130 - this is correct
# table(AllSampIndepMatrix)
# plot relatedness values
hist(AllSampIndepMatrix[AllSampIndepMatrix > 0 & AllSampIndepMatrix < 1],
     breaks = 100,
     xlab = "relatedness", 
     ylab = "pairs of individuals", 
     ylim = c(0, 2000),
     main = "Relatedness of independent individuals",
     col = "lightblue2")

##### RANDOM SAMPLING of 200 males and females
AllIndepFem <- subset(AllSamp, ID %in% AllSampIndep & sex == 2)
dim(AllIndepFem)
#head(AllIndepFem)

AllIndepMale <- subset(AllSamp, ID %in% AllSampIndep & sex == 1)
dim(AllIndepMale)
#head(AllIndepMale)

#set.seed(123)
#PoolMales <- sample(AllIndepMale$ID, size = 200, replace = FALSE, prob = NULL)
#length(PoolMales)
#PoolFemales <- sample(AllIndepFem$ID, size = 200, replace = FALSE, prob = NULL)
#length(PoolFemales)
PoolSeq <- read.csv("/Users/sonjalecic/Desktop/Sample_Lists/PoolSeqSamples.csv")


## check if there is overlap with the probelmatic (low, sheared DNA) samples
problematic <- read.csv("/Users/sonjalecic/Desktop/rereadyfortakeoff/ProblemDNASamples.csv")
intersect(PoolSeq$ID_males, problematic$ID) # no overlap with the problematic samples
intersect(PoolSeq$ID_females, problematic$ID) # no overlap with the probelmatic samples
