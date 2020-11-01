
##Fitness assay three Ace haplotypes - component fecundity
###egg count script
##Sonja_Lecic_06_2019


install.packages("lme4", "multcomp", "R.utils", "pastecs", "data.table")

library(lme4) 
library(multcomp)
library(R.utils)
library(pastecs)
library(data.table)

produce_counts_from_tables <- function(directory_path,
                                       file_extension=".xls", 
                                       user_bw=.2,
                                       min_px_Area=5,
                                       plot_dist_Area=FALSE,
                                       remove_small_outliers=TRUE,
                                       remove_large_outliers=TRUE,
                                       vec_names=NA){
  
  #user_bw: bandwidth used to compute the kernel density estimates of the Area distribution (see ?density)
  #min_px_Area: lower threshold that can be used to remove small dust particles. This threshold should not be set to the minimal expected Area of an egg (the first peak of a bimodal distribution is always removed in a later step)
  #plot_dist_Area: if TRUE, the function will plot the Area distribution with the cutoff used to remove outliers (Boolean)
  #remove_large_outliers: if TRUE, large particles are also removed from the final egg count
  #vec_names: each file name will be split using the "_" separator into several column to create a final data frame. The column names default values will be set to "ID_#" unless specified in vec_names.  
  
  temp_dir=getwd()
  setwd(directory_path)
  all_files=list.files(path=directory_path, pattern= file_extension)
  cat(sprintf("Loading input data from %s files... ", length(all_files)))
  read_all=lapply(all_files,read.table,sep="\t",header=TRUE)
  cat(sprintf("DONE\n"))
  
  counts=do.call(rbind,read_all)
  nb_lines_per_file=sapply(all_files,countLines)-1
  
  #Should be TRUE
  #sum(nb_lines_per_file)==nrow(counts)
  
  ## Add a column with the name of the file #id
  counts$id=rep(all_files, nb_lines_per_file)
  #counts=subset(counts,Area>= min_px_Area & counts$Major/counts$Minor < quantile(counts$Major/counts$Minor, .9999))
  counts=subset(counts,Area >= min_px_Area)
  
  density_dist=density(log10(counts$Area),bw=user_bw)
  tp=suppressWarnings(turnpoints(density_dist$y))
  
  # We use the first pit as the cutoff
  #(density_dist$x[tp$pits])[1]
  
  if(remove_small_outliers) {
    cat(sprintf("Remove %s small outliers (%s %%)\n", nrow(subset(counts,log10(Area)< density_dist$x[tp$peaks][1])),round(10000*nrow(subset(counts,log10(Area)< density_dist$x[tp$peaks][1]))/nrow(counts))/100))
  }
  
  cutoff_outliers=Inf
  if(remove_large_outliers){
    cutoff_outliers=density_dist$x[tp$peaks][2]+2*(density_dist$x[tp$peaks][2]-density_dist$x[tp$pits][1])
    cat(sprintf("Remove %s large outliers (%s %%)\n", nrow(subset(counts,Area> 10^cutoff_outliers)),round(1000*nrow(subset(counts,Area > 10^cutoff_outliers))/nrow(counts))/10))
  }
  
  ########
  
  if(plot_dist_Area){
    #quartz(height=3.5,width=6)
    plot(density_dist, xlab=expression("Area of the detected particles"~(log10(px^2))), ylab="Density", main="", bty = "n", lwd = 2)
    abline(v= (density_dist$x[tp$pits][1]),col="red")
    if(remove_large_outliers) abline(v= (cutoff_outliers),col="red")	
  }
  
  ## create the count table now
  counts_filt=subset(counts,log10(Area)< cutoff_outliers & log10(Area)>(density_dist$x[tp$pits][1]))
  
  labels_from_id=tstrsplit(strsplit(counts_filt$id,file_extension),"_")
  
  temp_data_frame=unique(data.frame(labels_from_id))
  temp_data_frame=cbind(temp_data_frame,tapply(counts_filt$Major, counts_filt$id, length))
  
  if(is.na(vec_names)){
    names(temp_data_frame)=c(paste0("ID_",1:length(labels_from_id)),"nb_eggs")
  }else{
    names(temp_data_frame)=c(vec_names,"nb_eggs")
  }
  
  setwd(temp_dir)
  
  return(list(temp_data_frame,counts_filt))
  
}

counts <- produce_counts_from_tables("/Users/sonjalecic/Desktop/Ace manuscript/Ace fitness assay/Fecundity/cold/")
counts[[1]] -> counts
counts %>% unite(cage, c(2:3), remove=FALSE) -> counts

### hot ###
hot = counts
hot_day1 = subset(hot, ID_4 == "Day 1")
hot_day1Iab = subset(hot_day1, !(ID_2 == "II"))
hot_Iab = subset(hot, !(ID_2 == "II"))

kruskal.test(hot_day1Iab$nb_eggs ~ hot_day1Iab$ID_2)
t.test(hot_day1Iab$nb_eggs ~ hot_day1Iab$ID_2)
wilcox.test(hot_day1Iab$nb_eggs ~ hot_day1Iab$ID_2)
boxplot(hot_day1$nb_eggs ~ hot_day1$ID_2, ylim = c(30, 700))
hist(hot_day1$nb_eggs)

model_0 = lmer(log(nb_eggs) ~ 1 + (1|ID_3), data=hot_Iab)
model_1 = lmer(log(nb_eggs) ~ ID_2 + (1|ID_3), data=hot_Iab)
anova(model_0, model_1)

hist(hot$nb_eggs)
model0 = glmer(nb_eggs ~ 1 + (1|ID_3) + (1|ID_4), data=hot_Iab, family = "poisson")

model1 = glmer(nb_eggs ~ ID_2 + ID_5 + (1|ID_3) + (1|ID_4), data=hot_Iab, family = "poisson")
anova(model0, model1)
summary(model1)
AIC(model1)

hot_cIab = subset(hot, !(ID_2 == "II"))
resaov = aov(nb_eggs ~ ID_2, data = hot_cIab)
summary(resaov)

#### hot ###
hot_Iab = subset(hot, !(ID_2 == "II"))
hist(hot$nb_eggs)
model0 = glmer(nb_eggs ~ 1 + (1|ID_3) + (1|ID_4), data=hot_Iab, family = "poisson")
model1 = glmer(nb_eggs ~ ID_2 + ID_5 + (1|ID_3) + (1|ID_4), data=hot_Iab, family = "poisson")
anova(model0, model1)
AIC(model1)

#### cold ###
cold_Iab = subset(cold, !(ID_2 == "II"))
hist(cold_Iab$nb_eggs)
model0 = glmer(nb_eggs ~ 1 + (1|ID_3) + (1|ID_4), data=cold_Iab, family = "poisson")
model1 = glmer(nb_eggs ~ ID_2 + ID_5 + (1|ID_3) + (1|ID_4), data=cold_Iab, family = "poisson")
anova(model0, model1)
AIC(model1)


## day 1 & 3 #
hot_day13 = subset(hot, !(ID_4 == "Day 2"))
boxplot(hot_day13$nb_eggs ~ hot_day13$ID_2, ylim = c(30, 700), col = c("grey", "red", "orange"))
boxplot(hot_day1$nb_eggs ~ hot_day1$ID_2, ylim = c(30, 700),  col = c("grey", "red", "orange"))

cold_am = subset(cold, ID_5 == "AM")
cold_amIab = subset(cold_am, !(ID_2 == "II"))
cold_Iab = subset(cold, !(ID_2 == "II"))
boxplot(cold_am$nb_eggs ~ cold_am$ID_2, ylim = c(30, 400), col = c("grey", "red", "orange"))

t.test(cold_Iab$nb_eggs ~ cold_Iab$ID_2)
kruskal.test(cold_Iab$nb_eggs ~ cold_Iab$ID_2)
chisq.test(cold_Iab$nb_eggs ~ cold_Iab$ID_2)



### cold ###
cold = counts
hist(cold$nb_eggs)
model0 = lmer(log(nb_eggs) ~ 1 + (1|ID_3) + (1|ID_4), data=cold)

model1 = lmer(log(nb_eggs) ~ ID_5 + (1|ID_3) + (1|ID_4), data=cold)
anova(model0, model1)


#### cold ####
cold_Ia = counts
cold_Ia_am = subset(cold_Ia, ID_5 == "AM")
boxplot(cold_Ia$nb_eggs, ylim = c(30, 400), main = "class Ia", col = "red")
boxplot(cold_Ia$nb_eggs ~ cold_Ia$ID_3, ylim = c(30, 400), main = "class Ia", col = "red")
boxplot(cold_Ia$nb_eggs ~ cold_Ia$ID_4, ylim = c(0, 400), main = "class Ia", col = "red")
boxplot(cold_Ia$nb_eggs ~ cold_Ia$ID_5, ylim = c(0, 400), main = "class Ia", col = "red")
# am #
cold_Ia_am = subset(cold_Ia, ID_5 == "AM")
boxplot(cold_Ia_am$nb_eggs, ylim = c(30, 400), main = "class Ia", col = "red")
# pm #
cold_Ia_pm = subset(cold_Ia, ID_5 == "PM")
boxplot(cold_Ia_pm$nb_eggs, ylim = c(30, 400), main = "class Ia", col = "red")
### day 2 ####
cold_Ia_day2 = subset(cold_Ia, ID_4 == "Day 2")
boxplot(cold_Ia_day2$nb_eggs, ylim = c(30, 400), main = "class Ia", col = "red")
boxplot(cold_Ia_day2$nb_eggs, cold_Ia$ID_5,  ylim = c(30, 400), main = "class Ia", col = "red")

cold_Ib = counts
boxplot(cold_Ib$nb_eggs, ylim = c(30, 400), main = "class Ia", col = "orange")
boxplot(cold_Ib$nb_eggs ~ cold_Ib$ID_3, ylim = c(30, 400), main = "class Ib", col = "orange")
boxplot(cold_Ib$nb_eggs ~ cold_Ib$ID_4, ylim = c(0, 400), main = "class Ib", col = "orange")
boxplot(cold_Ib$nb_eggs ~ cold_Ib$ID_5, ylim = c(0, 400), main = "class Ib", col = "orange")
# am #
cold_Ib_am = subset(cold_Ib, ID_5 == "AM")
boxplot(cold_Ib_am$nb_eggs, ylim = c(30, 400), main = "class Ib", col = "orange")
# pm #
cold_Ib_pm = subset(cold_Ib, ID_5 == "PM")
boxplot(cold_Ib_pm$nb_eggs, ylim = c(30, 400), main = "class II", col = "orange")
### day 2 ####
cold_Ib_day1 = subset(cold_Ib, ID_4 == "Day 1" )
boxplot(cold_Ib_day1$nb_eggs, ylim = c(30, 400), main = "class Ib", col = "orange")


#cold_II = counts
boxplot(cold_II$nb_eggs, ylim = c(30, 400), main = "class Ia", col = "grey")
boxplot(cold_II$nb_eggs ~ cold_II$ID_3, ylim = c(30, 400), main = "class II", col = "grey")
boxplot(cold_II$nb_eggs ~ cold_II$ID_4, ylim = c(0, 400), main = "class II", col = "grey")
boxplot(cold_II$nb_eggs ~ cold_II$ID_5, ylim = c(0, 400), main = "class II", col = "grey")
# am #
cold_II_am = subset(cold_II, ID_5 == "AM")
boxplot(cold_II_am$nb_eggs, ylim = c(30, 400), main = "class II", col = "grey")
# pm #
cold_II_pm = subset(cold_II, ID_5 == "PM")
# day 2 #
cold_II_day2 = subset(cold_II, ID_4 == "Day 2" )
boxplot(cold_II_day2$nb_eggs, ylim = c(30, 400), main = "class II", col = "grey")


## hot ##
hot_Ia = counts
hot_Ia = hot_Ia[-c(3,9), ]
add = data.frame(ID_1 = "hot", ID_2 = "Ia", ID_3 = "69", ID_4 = "Day 2", ID_5 = "AM", nb_eggs = 270)
ad2 = data.frame(ID_1 = "hot", ID_2 = "Ia", ID_3 = "13", ID_4 = "Day 2", ID_5 = "AM", nb_eggs = 270)
hot_Ia = rbind(hot_Ia, add, ad2)
boxplot(hot_Ia$nb_eggs, ylim = c(30, 700), main = "class Ia", col = "red")
boxplot(hot_Ia$nb_eggs ~ hot_Ia$ID_3, ylim = c(30, 700), main = "class Ia", col = "red")
boxplot(hot_Ia$nb_eggs ~ hot_Ia$ID_4, ylim = c(0, 700), main = "class Ia", col = "red")
boxplot(hot_Ia$nb_eggs ~ hot_Ia$ID_5, ylim = c(0, 700), main = "class Ia", col = "red")
# pm #
hot_Ia_pm = subset(hot_Ia, ID_5 == "PM")
boxplot(hot_Ia_pm$nb_eggs, ylim = c(30, 700), main = "class Ia", col = "red")
# am #
hot_Ia_am = subset(hot_Ia, ID_5 == "AM")
boxplot(hot_Ia_am$nb_eggs, ylim = c(30, 700), main = "class Ia", col = "red")
# day 2 #
hot_Ia_day2 = subset(hot_Ia, ID_4 == "Day 2" )
boxplot(hot_Ia_day2$nb_eggs, ylim = c(30, 400), main = "class Ia", col = "red")


hot_Ib = counts
boxplot(hot_Ib$nb_eggs, ylim = c(30, 700), main = "class Ib", col = "orange")
boxplot(hot_Ib$nb_eggs ~ hot_Ib$ID_3, ylim = c(30, 700), main = "class Ib", col = "orange")
boxplot(hot_Ib$nb_eggs ~ hot_Ib$ID_4, ylim = c(0, 700), main = "class Ib", col = "orange")
boxplot(hot_Ib$nb_eggs ~ hot_Ib$ID_5, ylim = c(0, 700), main = "class Ib", col = "orange")
# pm #
hot_Ib_pm = subset(hot_Ib, ID_5 == "PM")
boxplot(hot_Ib_pm$nb_eggs, ylim = c(30, 700), main = "class Ib", col = "orange")
# am #
hot_Ib_am = subset(hot_Ib, ID_5 == "AM")
boxplot(hot_Ib_am$nb_eggs, ylim = c(30, 700), main = "class Ib", col = "orange")
# day 2 #
hot_Ib_day2 = subset(hot_Ib, ID_4 == "Day 2")
boxplot(hot_Ib_day2$nb_eggs, ylim = c(30, 700), main = "class Ib", col = "orange")


hot_II = counts
boxplot(hot_II$nb_eggs, ylim = c(30, 700), main = "class II", col = "grey")
boxplot(hot_II$nb_eggs ~ hot_II$ID_3, ylim = c(30, 700), main = "class II", col = "grey")
boxplot(hot_II$nb_eggs ~ hot_II$ID_4, ylim = c(0, 700), main = "class II", col = "grey")
boxplot(hot_II$nb_eggs ~ hot_II$ID_5, ylim = c(0, 700), main = "class II", col = "grey")
# pm #
hot_II_pm = subset(hot_II, ID_5 == "PM")
boxplot(hot_II_pm$nb_eggs, ylim = c(30, 700), main = "class II", col = "grey")
# am #
hot_II_am = subset(hot_II, ID_5 == "AM")
boxplot(hot_II_am$nb_eggs, ylim = c(30, 700), main = "class II", col = "grey")





class_Ia_cold = data.frame(isoline = c(rep(13, 2), rep(69, 2)),
                           time = c(rep("morning", "evening", "morning", "evening")),
                           numeggs = c(51, 51, 155, 148))
boxplot(class_Ia_cold$numeggs, ylim = c(50, 400), main = "class Ia", col = "red")
boxplot(class_Ia_cold$numeggs, class_Ia_cold$isoline)

class_Ib_cold = data.frame(isoline = c(rep(26, 2), rep(75, 2)),
                           time = c(rep("morning", "evening", "morning", "evening")),
                           numeggs = c(84, 44, 338, 156))
boxplot(class_Ib_cold$numeggs, ylim = c(50, 400), main = "class Ib", col = "orange")
boxplot(class_Ib_cold$numeggs, class_Ib_cold$isoline)

class_II_cold = data.frame(isoline = c(rep(22, 4), rep(31, 4)),
                           time = c(rep("morning", "evening", "evening", "morning", "morning", "evening", "evening", "morning")),
                           numeggs = c(99, 187, 90, 68, 76, 52, 101, 37))
boxplot(class_II_cold$numeggs, ylim = c(50, 400), main = "class II", col = "grey40")
boxplot(class_II_cold$numeggs, class_II_cold$isoline)

par(mfrow = c(1, 3))


#### day 1 #####



##### cold morning + evening ####
class_Ia_cold = data.frame(isoline = c(rep(13, 3), rep(69, 3)),
                           time = c("morning", "morning", "evening", "morning", "morning", "evening"),
                           numeggs = c(51, 41, 51, 155, 61, 148))
boxplot(class_Ia_cold$numeggs, ylim = c(50, 400))
boxplot(class_Ia_cold$numeggs, class_Ia_cold$isoline)

class_Ib_cold = data.frame(isoline = c(rep(26, 3), rep(75, 3)),
                           time = c("morning", "morning", "evening", "morning", "morning", "evening"),
                           numeggs = c(84, 50, 44, 338, 116, 156))
boxplot(class_Ib_cold$numeggs, ylim = c(50, 400))
boxplot(class_Ib_cold$numeggs, class_Ia_cold$isoline)

class_II_cold = data.frame(isoline = c(rep(22, 3), rep(31, 3)),
                           time = c("morning", "morning", "evening", "morning", "morning", "evening"),
                           numeggs = c(99, 48, 184, 74, 59, 51))
boxplot(class_II_cold$numeggs, ylim = c(50, 400))
boxplot(class_II_cold$numeggs, class_Ia_cold$isoline)


#### hot morning + evening ####

class_Ia_hot = data.frame(isoline = c(rep(13, 3), rep(69, 3)),
                           time = c("morning", "morning", "evening", "morning", "morning", "evening"),
                           numeggs = c(621, 95, 501, 615, 64, 412))
boxplot(class_Ia_hot$numeggs, ylim = c(50, 700))
boxplot(class_Ia_hot$numeggs, class_Ia_cold$isoline)

class_Ib_hot = data.frame(isoline = c(rep(26, 3), rep(75, 3)),
                          time = c("morning", "morning", "evening", "morning", "morning", "evening"),
                          numeggs = c(261, 344, 146, 574, 437, 270))
boxplot(class_Ib_hot$numeggs, ylim = c(50, 700))
boxplot(class_Ib_hot$numeggs, class_Ia_cold$isoline)


class_II_hot = data.frame(isoline = c(rep(22, 3), rep(31, 3)),
                          time = c("morning", "morning", "evening", "morning", "morning", "evening"),
                          numeggs = c(387, 141, 501, 363, 330, 287))
boxplot(class_II_hot$numeggs, ylim = c(50, 700))
boxplot(class_II_hot$numeggs, class_Ia_cold$isoline)


#### hot morning + evening ####

class_Ia_hot = data.frame(isoline = c(rep(13, 4), rep(69, 4)),
                          time = c("morning", "evening", "evening", "morning", "morning", "evening", "evening", "morning"),
                          numeggs = c(643,  519, 209, 227, 628, 425, 249, 170))
boxplot(class_Ia_hot$numeggs, ylim = c(50, 700), main = "class Ia", col = "red")
boxplot(class_Ia_hot$numeggs, class_Ia_cold$isoline)

class_Ib_hot = data.frame(isoline = c(rep(26, 4), rep(75, 4)),
                          time = c("morning", "evening", "evening", "morning", "morning", "evening", "evening", "morning"),
                          numeggs = c(263,  146, 37, 194, 574,  270, 114, 275))
boxplot(class_Ib_hot$numeggs, ylim = c(50, 700), main = "class Ib", col = "orange")
boxplot(class_Ib_hot$numeggs, class_Ia_cold$isoline)


class_II_hot = data.frame(isoline = c(rep(22, 4), rep(31, 4)),
                          time = c("morning", "evening", "evening", "morning", "morning", "evening", "evening", "morning"),
                          numeggs = c(387,  500, 182, 250, 363, 287, 90, 256))
boxplot(class_II_hot$numeggs, ylim = c(50, 700), main = "class II", col = "grey40")
boxplot(class_II_hot$numeggs, class_Ia_cold$isoline)









