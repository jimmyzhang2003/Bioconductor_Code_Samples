library(ggplot2)

#load data
temp <- tempfile()
fileUrl <- "https://d396qusza40orc.cloudfront.net/rprog%2Fdata%2FProgAssignment3-data.zip"
download.file(fileUrl, temp)
hospitaldata <- read.csv(unz(temp, "outcome-of-care-measures.csv"), colClasses = "character")
unlink(temp)
dateDownloaded <- date()

#histogram of mortality rates
head(hospitaldata)
hospitaldata[, 11] <- as.numeric(hospitaldata[, 11])
mortalityrates <- hospitaldata[, 11]
g = ggplot(hospitaldata, aes(x = mortalityrates))
g = g + geom_histogram()
g = g + labs(x = "30-Day Mortality Rates from Heart Attack")
g = g + facet_wrap(.~ hospitaldata[, 7])
g


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
