library(fbi)
library(fredr)
library(tsutils)

subset <- snakemake@wildcards[["subset"]]
outfile <- snakemake@output[[1]]

####----- User Choices -----####

# Forecast horizon
h <- 1
# Degree of differening (dif = 2 adheres to Hauzenber et al., 2023)
dif <- 1
# Number of lags
lags <- 12

# Covariate choice
frednames <- c("INDPRO","IPMANSICS","CUMFNS","HOUST","PERMIT","USGOOD","MANEMP"
               ,"CES0600000007","AWHMAN","CUSR0000SAC","M1SL","NONBORRES","REALLN",
               "COMPAPFFx","TB3SMFFM","TB6SMFFM","T1YFFM","T10YFFM","AAAFFM","BAAFFM")

####----- Preparation Steps -----####

## 1: Supply FRED mnemonics for data download.  We follow Hauzenber et al. (2023)
filepath <- "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2023-03.csv"

startdate <- as.Date("01/01/1980","%m/%d/%Y")
enddate <- as.Date("01/01/2023","%m/%d/%Y")

data <- fredmd(filepath, date_start = startdate, date_end = enddate, transform = TRUE)

## 2: Format data into a data frame
datevec <- data$date

if (subset == "subset") {
  data <- data[,frednames]
} else {
  data <- data[,2:ncol(data)]
}

## 3: Impute any missing values following Ercument Cahan, Jushan Bai, and Serena Ng (2021)
data_est <-  tp_apc(data, 4, center = FALSE, standardize = FALSE, re_estimate = TRUE)
data_imp <- data_est$data

## 3: Get inflation data
fredr_set_key("315410e54f1b6f2552a99cefd47e2344") #API key
inflation <- fredr(  series_id = "CPIAUCSL",
                     observation_start = startdate,
                     observation_end = enddate,
                     frequency = "m") #percent change transformation 
inflation <- inflation$value

if (dif == 2) {
  inflation <- log(inflation[(h+2):length(inflation)]/inflation[2:(length(inflation)-h)]) - log(inflation[(2):(length(inflation)-h)]/inflation[1:(length(inflation)-h-1)])
} else if (dif ==1) {
  inflation <- log(inflation[(h+2):length(inflation)]/inflation[2:(length(inflation)-h)])*100 # following Chan (2017) The Stochastic Volatility in Mean Model With TimeVarying Parameters
}

X <- as.matrix(data_imp[2:(dim(data_imp)[1]-h),])
T <- length(inflation)
if (subset == "subset") {
  K <- length(frednames)
} else {
  K <- ncol(X)
}
# Create matrix of lags
X_all <- array(0,c(T-lags,K*lags))

for (j in 1:K){
  lagtemp <- lagmatrix(X[,j],1:lags)
  X_all[,(lags*(j-1)+1):(lags*j)] <- lagtemp[((lags+1):dim(lagtemp)[1]),]
}

####----- Save Data -----####

# placeholder for data matrix
y <- as.vector(inflation[(lags+1):length(inflation)])
datevec <- datevec[(lags+1):length(inflation)]
if (subset == "subset") {
  Xnames <- frednames
} else {
  Xnames <- colnames(data_imp[2:ncol(data_imp)])
}
yname <- "CPIAUCSL"
lagstructure <- rep(1:lags,K)

emp_app <- list(X = X_all,
                y = y,
                Xnames = Xnames,
                yname = yname,
                datevec = datevec,
                lagstructure = lagstructure)

saveRDS(emp_app, file = outfile)
