library(ggplot2)

#load data
temp <- tempfile()
fileUrl <- "https://d396qusza40orc.cloudfront.net/rprog%2Fdata%2FProgAssignment3-data.zip"
download.file(fileUrl, temp)
hospitaldata <- read.csv(unz(temp, "outcome-of-care-measures.csv"), colClasses = "character")
unlink(temp)
dateDownloaded <- date()

#histogram of heart attack mortality rates
head(hospitaldata)
hospitaldata[, 11] = as.numeric(hospitaldata[, 11])
hattackrate = hospitaldata[, 11]
h = ggplot(hospitaldata, aes(x = hattackrate))
h = h + geom_histogram()
h = h + labs(x = "30-Day Mortality Rates from Heart Attack")
h = h + facet_wrap(.~ hospitaldata[, 7])
h

#boxplot of heart failure mortality rates
hospitaldata[, 17] = as.numeric(hospitaldata[, 17])
hfailrate = hospitaldata[, 17]
b = ggplot(hospitaldata, aes(x = hospitaldata[, 7], y = hfailrate))
b = b + geom_boxplot()
b = b + labs(x = "State", 
             y = "30-Day Mortality Rates from Heart Failure")
b

#violin plot of pneumonia mortality rates
hospitaldata[, 23] = as.numeric(hospitaldata[, 23])
pneumorate = hospitaldata[, 23]
p = ggplot(hospitaldata, aes(x = hospitaldata[, 7], y = pneumorate))
p = p + geom_violin()
p = p + labs(x = "State",
             y = "30-Day Mortality Rates from Pneumonia")
p

#function for best hospital in a state
best <- function(state, outcome) {
        if(!exists("hospitaldata")) {
          temp <- tempfile()
          fileUrl <- "https://d396qusza40orc.cloudfront.net/rprog%2Fdata%2FProgAssignment3-data.zip"
          download.file(fileUrl, temp)
          hospitaldata <<- read.csv(unz(temp, "outcome-of-care-measures.csv"), colClasses = "character")
          unlink(temp)  
        }
        if(state %in% hospitaldata[, 7] == FALSE) {
            stop("invalid state")
        } else if(outcome == "heart attack") {
            rows <- which(hospitaldata[, 7] == state)
            df <- data.frame(hospitaldata[rows, c(2, 11)])
            bestrow <- which(grepl(min(df[, 2]), df[, 2]))
            besthospital <- df[bestrow, 1]
            print(besthospital)
        } else if(outcome == "heart failure") {
            rows <- which(hospitaldata[, 7] == state)
            df <- data.frame(hospitaldata[rows, c(2, 17)])
            bestrow <- which(grepl(min(df[, 2]), df[, 2]))
            besthospital <- df[bestrow, 1]
            print(besthospital)
        } else if(outcome == "pneumonia") {
            rows <- which(hospitaldata[, 7] == state)
            df <- data.frame(hospitaldata[rows, c(2, 23)])
            bestrow <- which(grepl(min(df[, 2]), df[, 2]))
            besthospital <- df[bestrow, 1]
            print(besthospital)
        } else {
            stop("invalid outcome") 
        }
      }
