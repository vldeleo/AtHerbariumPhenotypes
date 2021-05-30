
# I don't think I can 100% guarantee Zeros are correct, although I have checked some of them and fixed errors where negative area was calculated

#Photothermal Units
library(pracma)
library(geosphere)
library(rootSolve)

# cleaning to get rid of accessions that have missing latitude or temperature data
FTclim <- read.csv('fileaddress/PhenoLit_Aug27_clim.csv', stringsAsFactors = F) #load climate data
FTlit <- read.csv("fileaddress/PhenoLit_Aug26_ele.csv", stringsAsFactors = FALSE) #load data from literature that includes flowering date
FTclim$JulDay <- FTlit$JulDay[match(FTclim$Ecotype, FTlit$Ecotype)] #need date of collection
for (c in 5:40){
  FTclim[,c] <- FTclim[,c]/10 #for whatever reason, worldclim temps are reported in tenths of degrees, so need to divide
}
for (i in 1:nrow(FTclim)){
  FTclim$index[i] <- i
} # and index to ensure I can match my ptu to my accessions despite any errors in name

use <- subset(FTclim, !is.na(FTclim$Latitude) & FTclim$Latitude <= 90) 
lousy <- subset(use, is.na(use$tmn_1) & is.na(use$tmn_2) & is.na(use$tmn_3)) #too much missing climate data makes the calculation impossible
almost <- use[!(rownames(use) %in% rownames(lousy)),]
good <- almost[(almost$JulDay != 1),] #these collection dates are not accurate
good <- good[which(good$Latitude < 65.72),]
#^ this removes the accessions where there isn't a dawn or dusk in summer months. these ecotypes can be calculated separately in the highlat script.


dayofmean <- c(15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319,349.5,380.5)
#a very rough assumption that the mean temperature of the month should fall on middle day of the month.

#function to find hourly temperature in a day, using Burghardt et al. equation
##maybe don't need daymin/daymax as separate functions since I only need them in this loop
degreeday <- function(d, x) {
  result <- vector('numeric', 24L)
  daymin <- function(x, d, dayofmean) {
    inside <- splinefun(x = dayofmean, y = good[x, c(which(names(good) == "tmin_1"):which(names(good) == "tmin_12"),which(names(good) == "tmin_1"))], method = 'periodic') 
    return(inside(d))}
  daymax <- function(x, d, dayofmean) {
    inside = splinefun(x = dayofmean, y = good[x, c(which(names(good) == "tmax_1"):which(names(good) == "tmax_12"),which(names(good) == "tmax_1"))], method = 'periodic')
    return(inside(d))}
  Hn <- (24-daylength(good[x,'Latitude'], d))/2
  Hs <- Hn + daylength(good[x,'Latitude'], d)
  Tn <- daymin(x, d, dayofmean)
  Tm <- daymax(x, d, dayofmean)
  Tp <- daymin(x, d+1, dayofmean)
  Ts <- Tm - .227538 * (Tm - Tp)
  Hp <- ((24-daylength(good[x,'Latitude'], d+1))/2)
  Hsb <- 24 - (24-daylength(good[x,'Latitude'], d-1))/2 
  Hm <- (Hn + 2.036391 * sin(2 * pi * (d-79.22015)/365) + 9.285504)
  
  for (t in 1:24){
    result[t] <- ifelse( Hm >= t & t >= Hn, (Tm + Tn)/2 + (Tm -Tn)/2 * cos(pi + (t - Hn)/(Hm-Hn) * pi), ifelse( Hs >= t & t > Hm, (Ts + (Tm-Ts) * log(((1 + Hs - Hm)-(t - Hm)), base = (1 + Hs - Hm))), ifelse( t > Hs, ((Tp - Ts)/sqrt(Hp + 24 - Hs) * sqrt((t-Hs)) + Ts), ((Tp - Ts)/sqrt(abs(Hp + 24 - Hs)) * sqrt((t + 24 - Hsb)) + Ts))))
    }
  return(data.frame(result))
}

#function to extend hourly temperatures calculations to the whole year
dhot <- function(x) {
  res <- vector('numeric')
  for (d in 1:365){
    res[d] <- degreeday(d, x)
  }
  temp <- stack(as.data.frame(res))
  object <- c(temp$values)
  return(object)
}

final <- vector()
solution <- vector()
solution30 <- vector() 
names <- vector()
orivalue <- vector()
oriname <- vector()

#in case my program gets interrupted partway through:
#solution <- read.csv('fileaddress/...')
#solution <- solution[,1]


#start of PTU calculations
for (x in 1:nrow(good)){
  message(paste('Iteration ',x,' of ',nrow(good), sep = ''))
  Zeros30a <- vector()
  Zeros30b <- vector()
  area <- vector()
  area30 <- vector()
  #these need to be cleared at each loop, or accessions that are collected earlier than accessions that come before them in the dataset will be messed up
  
  #each hour of daybreak and sunset until day of collection  
  dawn <- vector()
  dawntally <- vector()
  dusk <- vector()
  dusktally <- vector()
  for (d in 1:good[x, 'JulDay']) {
    dawn[d] <- (12 - daylength(good[x,'Latitude'], d)/2)
    dawntally[d] <- (dawn[d] + (d-1)*24)
  }
  for (d in 1:good[x, 'JulDay']) {
    dusk[d] <- (12 + daylength(good[x,'Latitude'], d)/2)
    dusktally[d] <- (dusk[d] + (d-1)*24)
  }
  
  #find every temperature at every hour
  p <- dhot(x)
  
  #turn the hourly temperatures into a curve that I can calculate area under
  Equa <- splinefun(x = 1:8760, y = p, method = 'natural')
  
  equamin <- function(h) {
    Res <- (Equa(h) - 3) #I set the minimum temperature of growth to 3 degrees, because At isn't thriving at low temps. the biological minimum might be higher than this, or even accession dependent.
    return(Res)
  }
  
  end <- (good[x, 'JulDay'] * 24)  
  
  # do I have to find zeros for this accession?
  peak <- vector()
  trough <- vector()
  for (h in 1:good[x, 'JulDay']){
    peak[h] <- max(equamin(dawntally[h]:dusktally[h]))
    trough[h] <- min(equamin(dawntally[h]:dusktally[h]))
  }
  
  w <- (peak > 0 & trough < 0)
  
  if (sum(w) == 0) {
    if (max(equamin(0:end)) < 0) { 
      solution[x] <- NA
    } else {
      for (i in 1:good[x, 'JulDay']){
        area[i] <- quad(equamin, min(dawntally[i]), max(dusktally[i]))
        solution[x] <- sum(area)
      }
    }
  } else {
    
    ##create vector of midnights
    daydemarc <- vector()
    for (d in 1:(good[x, 'JulDay']+1)){
      daydemarc[d] <- ((d-1) *24)
    }  
    ##create extra dawn in case Zeros2 is greater than last dusk
    dawntally[(good[x, 'JulDay']+1)] <- (12 - (daylength(good[x,'Latitude'], (good[x, 'JulDay']+1))/2) + end)
    
    
    #have to find hours where temperature climbs above 3 degrees to compare back to dawn and dusk
    ZerosAll <- vector()
    for (i in 1:(length(daydemarc)-1)) {
      sol <- uniroot.all(equamin, lower = daydemarc[i], upper = daydemarc[i+1], tol = .000001)
      ZerosAll <- c(ZerosAll, sol)
    }
    # just in case the temperature calculation starts above 0 and drops
    #if(ZerosAll[1] < dawntally[1]){
    #  ZerosAll <- ZerosAll[-1]
    #}
    # note, this would work if the equation ALWAYS makes dawn the coldest/near coldest
    #   time of day, so there will not be a 0 before dawn unless temperature is falling
    #   but I could be wrong about this

    
    Zeros <- vector()
    Zeros2 <- vector()
    
    q = min(which(trough < 0 & peak > 0))
      if((trough[q] > trough[q+1]) & (equamin(0) > 0)){  #if the temperature is decreasing, and the temperature starts above 0, the first Zero will need to be treated as a potential endpoint
      Zeros2 <- ZerosAll[seq(1, length(ZerosAll), by = 2)]
      Zeros <- ZerosAll[seq(2, length(ZerosAll), by = 2)]
    } else {
      if((trough[q] < trough[q+1]) & (equamin(0) > 0)){ #if the temperature is increasing, and starts above 0, but still dips below 0 at some point, the first Zero will need to be treated as a potential endpoint
        Zeros2 <- ZerosAll[seq(1, length(ZerosAll), by = 2)]
        Zeros <- ZerosAll[seq(2, length(ZerosAll), by = 2)]
        #in hindsight, perhaps both the above cases can be captured more simply as if(equamin(0) > 0)
        } else{
        Zeros2 <- ZerosAll[seq(2, length(ZerosAll), by = 2)]
        Zeros <- ZerosAll[seq(1, length(ZerosAll), by = 2)]  
      }
    }
  
    
      
    # find Zeros and dawn values in each day and compare
    # want greater value for PTU calculation, but don't want values for days < 3
    ori <- vector()
    for (i in 1:(length(daydemarc)-1)) {
      if (sum((Zeros > daydemarc[i]) & (Zeros < daydemarc[i +1]))  > 0)  { #if the temperature crosses 0 at some point during the day
        if (dawntally[which((dawntally > daydemarc[i]) & (dawntally < daydemarc[i +1]))] > Zeros[which((Zeros > daydemarc[i]) & (Zeros < daydemarc[i + 1]))]) { #if the temperature is above 3* by sun-up
          ori[i] <- dawntally[which((dawntally > daydemarc[i]) & (dawntally < daydemarc[i +1]))] #use dawn as the starting time for calculation
        } else {
          ori[i] <- Zeros[which((Zeros > daydemarc[i]) & (Zeros < daydemarc[i + 1]))] #if the temperature is not above 3* by sun-up, use the temperature threshold as the starting time
          if (Zeros[(Zeros > daydemarc[i]) & (Zeros < daydemarc[i +1])] > dusktally[i]){  #unless the temperature gets above 3* AFTER sundown, then don't calculate ANY PTUs (no start point for integration)
            ori[i] <- NA 
          }
          }
      } else {
        ori[i] <- dawntally[which((dawntally > daydemarc[i]) & (dawntally < daydemarc[i +1]))] #if the temperature doesn't cross zero all day, use dawn as starting time
        if (max(equamin(dawntally[i]:dusktally[i])) < 0 ){ #(unless it doesn't cross 0 because it's too low! then, start time is NA)
          ori[i] <- NA
        }}
    }
    
    
    #for some accessions, the temperature gets above 0 for such a tiny amount of time one day that the R calculation to find zeros doesn't find two intercepts in the interval
    #because this is such a tiny positive area, and because I don't want the negative value to skew my result, I can leave off the first ori point in these. don't use this option without visualizing the PTU curve first.
    
#    if(x %in% c(79, 463, 464)){
#      orivals <- which(!is.na(ori))
#      badori <- min(orivals)
      
#      ori[badori] <- NA
 #     }
    #I wanted to replace the above with a more inclusive line 194: max(equamin(dawntally[i]:dusktally[i])) < 0 | (max(equamin(dawntally[i]:dusktally[i])) > 0 & min(equamin(dawntally[i]:dusktally[i])) < 0)
      #but this messed up accessions where equamin could find no 0s, but the temperature in fact briefly dropped a tiny bit below 0. 
    
    
    # find nearest Zeros2 and nearest dusktally values and compare
    # want lesser value for PTU calculation, but don't want values for days < 3 degrees
    term <- vector()
    for (i in 1:(length(daydemarc)-1)) {
      #I use dawntally instead of daydemarc here because the temperature decreases up until dawn and I don't want to include Zeros2 that happen before dawn
      if (sum((Zeros2 > dawntally[i]) & (Zeros2 < dawntally[i+1])) == 1){ #if the temperature crosses 3* at some point before dawn the next day
        if (dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))] > Zeros2[which((Zeros2 > dawntally[i]) & (Zeros2 < dawntally[i+1]))]) { #and if dusk is later than when the temperature dips down (but the temperature dip is before dawn)
          term[i] <- Zeros2[which((Zeros2 > dawntally[i]) & (Zeros2 < dawntally[i+1]))] #then use zero2 as the endpoint for calculationg PTUs
          } else {
          term[i] <- dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))]
          if(sum(Zeros[which((Zeros > daydemarc[i]) & (Zeros < daydemarc[i +1]))]) > 0 && Zeros[which((Zeros > daydemarc[i]) & (Zeros < daydemarc[i +1]))] > dusktally[i] &  Zeros2[which((Zeros2 > dawntally[i]) & (Zeros2 < dawntally[i+1]))] > dusktally[i]) { #unless both zeros fall after dusk, then don't calculate and endpoint
            term[i] <- NA
          }
          }  
        if (!is.na(ori[i]) & Zeros2[which((Zeros2 > dawntally[i]) & (Zeros2 < dawntally[i+1]))] < ori[i]){
          term[i] <- dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))]
        }
        # I don't really know how ^ would be possible after deleting extra ZerosAll
        if (sum((Zeros > daydemarc[i]) & (Zeros < daydemarc[i+1])) < 1){
          term[i] <- dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))]
        }
      } else { 
        term[i] <- dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))]
        if (is.na(ori[i])) {
          term[i] <- NA
        }
      }
    }
    
    ori2 <- ori[!is.na(ori)]
    term2 <- term[!is.na(term)]
    
    for (i in 1:length(term2)){
      area[i] <- quad(equamin, min(ori2[i]), max(term2[i]))
    }
    if (any(area < 0)) message("Error in Zero Assignment")
        solution[x] <- sum(area)
  
}  
  
  ## now I calculate time above 30 degrees and delete extra area
  ##### I did not use the sub-30 PTU calculation in my paper, so I am less confident it is error-free   ####
  ## I didn't find a big difference between the two calculations

  if (max(Equa(0:end)) > 30){
    message('higher than 30')
    equamax <- function(h){
      Res <- (Equa(h) - 30)
      return(Res)  
    }
    
    daydemarc <- vector()
    for (d in 1:(good[x, 'JulDay']+1)){
      daydemarc[d] <- ((d-1) *24)
    } 
    
    #need extra dawn for calculating ori 
    #I did this up in the sub30 calculation, do I need to repeat it?
    for (d in 1:(good[x, 'JulDay']+1)) {
      dawn[d] <- 12 - daylength(good[x,'Latitude'], d)/2
      dawntally[d] <- (dawn[d] + (d-1)*24)
    }
    
    ZerosAll30 <- vector()
    for (i in 1:(length(daydemarc)-1)) {
      sol <- uniroot.all(equamax, lower = daydemarc[i], upper = daydemarc[i+1], tol = .0001)
      ZerosAll30 <- c(ZerosAll30, sol)
    }

      
    Zeros30a <- vector()
    Zeros30b <- vector()
    
    peak30 <- vector()
    trough30 <- vector()
    for (h in 1:good[x, 'JulDay']){
      peak30[h] <- max(equamax(dawntally[h]:dusktally[h]))
      trough30[h] <- min(equamax(dawntally[h]:dusktally[h]))
    }
    q30 <- min(which(trough30 < 0 & peak30 > 0))
    
    if(trough30[q30] > trough30[q30+1] & equamax(0) > 0){
      Zeros30b <- ZerosAll30[seq(1, length(ZerosAll30), by = 2)]
      Zeros30a <- ZerosAll30[seq(2, length(ZerosAll30), by = 2)]
    } else {
      if((trough30[q30] < trough30[q30+1]) & (equamax(0) > 0)){
        Zeros30b <- ZerosAll30[seq(1, length(ZerosAll30), by = 2)]
        Zeros30a <- ZerosAll30[seq(2, length(ZerosAll30), by = 2)]
      } else{
        Zeros30b <- ZerosAll30[seq(2, length(ZerosAll30), by = 2)]
        Zeros30a <- ZerosAll30[seq(1, length(ZerosAll30), by = 2)]  
      }
      }
    
    
    ori30 <- vector()
    for (i in 1:(length(daydemarc)-1)) {
      if (sum((Zeros30a > daydemarc[i]) & (Zeros30a < daydemarc[i +1])) > 0)  {
        if (dawntally[which((dawntally > daydemarc[i]) & (dawntally < daydemarc[i +1]))] > Zeros30a[which((Zeros30a > daydemarc[i]) & (Zeros30a < daydemarc[i + 1]))]) {
          ori30[i] <- dawntally[which((dawntally > daydemarc[i]) & (dawntally < daydemarc[i +1]))]
        } else {
          ori30[i] <- Zeros30a[which((Zeros30a > daydemarc[i]) & (Zeros30a < daydemarc[i + 1]))]
          if (Zeros30a[(Zeros30a > daydemarc[i]) & (Zeros30a < daydemarc[i +1])] > dusktally[i]){  #unless the temperature gets above 30* AFTER sundown, then don't calculate ANY PTUs (no start point for integration)
            ori30[i] <- NA 
          }
        }
      } else {
        ori30[i] <- dawntally[which((dawntally > daydemarc[i]) & (dawntally < daydemarc[i +1]))]
        if ((Zeros30a[1] > daydemarc[i]) & !(Zeros30a[1] < daydemarc[i +1]) & (Zeros30b[1] > Zeros30a[1])) {
          ori30[i] <- NA
        }}
    }
    
    
    term30 <- vector()
    for (i in 1:(length(daydemarc)-1)) {
      if (sum((Zeros30b > dawntally[i]) & (Zeros30b < dawntally[i+1])) == 1){
        if (dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))] > Zeros30b[which((Zeros30b > dawntally[i]) & (Zeros30b < dawntally[i+1]))]) {
          term30[i] <- Zeros30b[which((Zeros30b > dawntally[i]) & (Zeros30b < dawntally[i+1]))]
        } else {
          term30[i] <- dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))]
          
        if (Zeros30b[which((Zeros30b > dawntally[i]) & (Zeros30b < dawntally[i+1]))] < ori30[i]){
          term30[i] <- dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))]
        }
        if (sum((Zeros30b > daydemarc[i]) & (Zeros30b < daydemarc[i+1])) < 1){
          term30[i] <- dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))]
        }
        }
      } else { 
        term30[i] <- dusktally[which((dusktally > daydemarc[i]) & (dusktally < daydemarc[i +1]))]
        if (is.na(ori30[i])) {
          term30[i] <- NA
        }
      }
    }
    
    ori30 <- ori30[!is.na(ori30)]
    term30 <- term30[!is.na(term30)]
    
    for (i in 1:length(term30)){
      area30[i] <- quad(equamax, min(ori30[i]), max(term30[i]))
    }
    if (any(area30 < 0)) message("Error in Zero Assignment > 30")
    solution30[x] <- sum(area30)
    final[x] <- (solution[x] - solution30[x])

     } else {
    final[x] <- solution[x]
     }
  
#I find it helpful to visually check that a valid start/stop point was calculated for each accession
  #if you are not running this code interactively, skip the block below
  plot(equamin, (min(ori2, na.rm = T) - 50), (min(ori2, na.rm = T) + 250), n = 500, main = c("Accession ", x))
  abline(h=0)
  abline(v= dawntally, col = 'gold')
  abline(v= dusktally, col = 'navy blue')
  points(term, equamin(term), col = 'red')
  points(ori, equamin(ori), col = 'green') #plot the term point first, so it doesn't hide whether there was a valid origin for integration
  
  orivalue[x] <- as.character(good$index[x])
  oriname[x] <- as.character(good$Ecotype[x])
  solution <- round(solution, digits = 4)
  final <- round(final, digits = 4)
  
  final2 <- cbind.data.frame(orivalue, oriname, final, stringsAsFactors = FALSE)
  solution2 <- cbind.data.frame(orivalue, oriname, solution, stringsAsFactors = FALSE)
  ### end of >30 calculations ###
  ###make sure to set final in the dataframes to be numeric
  
  
  write.csv(solution2, 'fileaddres/PTU.csv', row.names = FALSE)  
  write.csv(final2, 'fileaddress/PTUsub30.csv', row.names = FALSE)
}
#}


  