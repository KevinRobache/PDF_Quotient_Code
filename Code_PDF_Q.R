# This function computes the 1D probability density function (PDF) of a time series x.
# The bandwidth is defined as a fraction of the data range 
# (e.g., range_bw = 0.01 corresponds to 1 % of the range between min(x) and max(x)).
# It returns a DataFrame containing the counts and probabilities for each interval.
PDFfun <- function(x, range_bw){
  min_x <- min(x, na.rm = T)
  max_x <- max(x, na.rm = T)
  bw <- (max_x - min_x) * range_bw
  seq_x <- seq(min_x, max_x, bw)
  for (i in 2:length(seq_x)){
    if (i == 2){
      Count <- length(x[which(x >= seq_x[i-1] & x <= seq_x[i])])
      Range <- (seq_x[i-1] + seq_x[i])/2
    }
    else{
      Count <- c(Count, length(x[which(x > seq_x[i-1] & x <= seq_x[i])]))
      Range <- c(Range, (seq_x[i-1] + seq_x[i])/2)
    }
  }
  Proba <- Count/sum(Count)
  Results_PDF <- data.frame(Range, Count, Proba)
  return(Results_PDF)
}


# This function computes the 2D probability density function (PDF) of two time series, x and y.
# The bandwidth is defined as a fraction of the data ranges:
#   - for x: range_bw * (max(x) - min(x))
#   - for y: range_bw * (max(y) - min(y))
# For example, range_bw = 0.01 corresponds to 1 % of the respective ranges of x and y.
# The function returns a DataFrame containing the counts and probabilities for each (x, y) interval.
PDF_2D_fun <- function(x, y, range_bw){
  min_x <- min(x, na.rm = T); max_x <- max(x, na.rm = T)
  min_y <- min(y, na.rm = T); max_y <- max(y, na.rm = T)
  bw_x <- (max_x - min_x) * range_bw; bw_y <- (max_y - min_y) * range_bw
  seq_x <- seq(min_x, max_x, bw_x); seq_y <- seq(min_y, max_y, bw_y)
  data <- data.frame(x,y)
  for (i in 2:length(seq_x)){
    for (j in 2:length(seq_y)){
      if (i == 2 & j == 2){
        Count <- nrow(data[which((data$x >= seq_x[i-1] & data$x <= seq_x[i]) & (data$y >= seq_y[j-1] & 
                                                                                  data$y <= seq_y[j])),])
        Range_x <- (seq_x[i-1] + seq_x[i])/2
        Range_y <- (seq_y[j-1] + seq_y[j])/2
      }
      else{
        Count <- c(Count, nrow(data[which((data$x >= seq_x[i-1] & data$x <= seq_x[i]) & 
                                            (data$y >= seq_y[j-1] & data$y <= seq_y[j])),]))
        Range_x <- c(Range_x, (seq_x[i-1] + seq_x[i])/2)
        Range_y <- c(Range_y, (seq_y[j-1] + seq_y[j])/2)
      }
    }
  }
  Proba <- Count/sum(Count)
  Results_PDF <- data.frame(Range_x, Range_y, Count, Proba)
  return(Results_PDF)
}


# This function computes the PDF quotient of two time series, x and y, 
# using the previously defined 1D and 2D PDF functions.
#
# The bandwidth is defined as a fraction of the data ranges:
#   - for x: range_bw * (max(x) - min(x))
#   - for y: range_bw * (max(y) - min(y))
# For example, range_bw = 0.01 corresponds to 1% of the respective ranges of x and y.
#
# The function returns a DataFrame containing the counts, probabilities, 
# and the quotient Q for each (x, y) interval.
Q_ratio <- function(x, y, range_bw){
  Ind <- which(!is.na(x) & !is.na(y)); x <- x[Ind]; y <- y[Ind]
  PDF_x <- PDFfun(x, range_bw = range_bw); PDF_y <- PDFfun(y, range_bw = range_bw)
  PDF_2D <- PDF_2D_fun(x, y, range_bw = range_bw)
  PDF_2D <- subset(PDF_2D, Proba != 0)
  PDF_x <- subset(PDF_x, Proba != 0); PDF_x <- subset(PDF_x, Proba != 0)
  colnames(PDF_x)[1] <- 'Range_x'; colnames(PDF_y)[1] <- 'Range_y'
  colnames(PDF_2D)[c(3,4)] <- c('Count_2D', 'Proba_2D')
  colnames(PDF_x)[c(2,3)] <- c('Count_x', 'Proba_x')
  colnames(PDF_y)[c(2,3)] <- c('Count_y', 'Proba_y')
  PDF_final <- dplyr::full_join(PDF_2D, PDF_x, by = 'Range_x')
  PDF_final <- dplyr::full_join(PDF_final, PDF_y, by = 'Range_y')
  PDF_final <- PDF_final[order(PDF_final$Range_x, decreasing = F),]
  PDF_final$Q <- PDF_final$Proba_2D / (PDF_final$Proba_x * PDF_final$Proba_y)
  PDF_final$Q <- log10(PDF_final$Q); PDF_final <- subset(PDF_final, !is.infinite(Q))
  return(PDF_final)
}
