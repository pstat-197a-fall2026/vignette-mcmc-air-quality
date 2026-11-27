###############################################################################
#  01_preprocess.R
#  Combine the 3 PRSA station files into one clean, model-ready daily dataset.
#
#  Input : data/PRSA_Data_<Station>_20130301-20170228.csv  (3 files, hourly)
#  Output: data/air_daily.csv   (one row per station per day)
#          data/scaling.csv     (means/SDs used to standardize predictors)
###############################################################################

data_dir <- "/Users/jocelynlee/Downloads/vignette-mcmc-air-quality/data"
files <- list.files(data_dir, pattern = "PRSA_Data_.*\\.csv$", full.names = TRUE)
stopifnot(length(files) == 3)

# ---- Step 1: read and stack ------------------------------------------------
raw <- do.call(rbind, lapply(files, read.csv))
cat("Raw rows:", nrow(raw), "\n")

# ---- Step 2: keep only the columns we need ---------------------------------
#   PM2.5 = response (ug/m^3), WSPM = wind speed (m/s),
#   RAIN  = precipitation (mm),  PRES = air pressure (hPa)
raw <- raw[, c("station", "year", "month", "day", "hour",
               "PM2.5", "WSPM", "RAIN", "PRES")]

# James I do not know if you want to impute but then we could change this line 
# # ---- Step 3: check missing values ------------------------------------------
# cat("\nMissing values per column:\n")
# print(colSums(is.na(raw)))
# # A small share of hours are missing. We drop incomplete hours rather than
# # impute; with ~35,000 hours per station this loses very little.
# raw <- raw[complete.cases(raw), ]
# cat("Rows after dropping incomplete hours:", nrow(raw), "\n")

# ---- Step 4: aggregate hourly -> daily -------------------------------------
# Why: consecutive hours are almost copies of each other (strong
# autocorrelation). Regression treats every row as independent information,
# so hourly data would make our uncertainty look far smaller than it is.
raw$date <- as.Date(sprintf("%d-%02d-%02d", raw$year, raw$month, raw$day))

daily <- aggregate(cbind(PM2.5, WSPM, PRES) ~ station + date, data = raw, FUN = mean)
rain  <- aggregate(RAIN ~ station + date, data = raw, FUN = sum)   # daily TOTAL rain
n_hrs <- aggregate(hour ~ station + date, data = raw, FUN = length)
names(n_hrs)[3] <- "n_hours"

daily <- merge(merge(daily, rain), n_hrs)
daily <- daily[daily$n_hours >= 18, ]     # drop days missing > 6 hours of data

# ---- Step 5: transform the response ----------------------------------------
# PM2.5 is strongly right-skewed (a few extreme smog days), which violates the
# normal-errors assumption for linear regression. log() makes it roughly symmetric. Coefficients then
# read as approximate % changes: exp(beta) - 1.
daily <- daily[daily$PM2.5 > 0, ]
daily$log_pm25 <- log(daily$PM2.5)

# ---- Step 6: standardize predictors (mean 0, SD 1) -------------------------
# (a) Coefficients become comparable: "effect of a 1-SD change".
# (b) HMC samples much better when parameters are on similar scales.
#     (Pressure is ~1000 hPa; unscaled it would wreck the sampler's step size.)
preds <- c("WSPM", "RAIN", "PRES")
scaling <- data.frame(variable = preds,
                      mean = sapply(daily[preds], mean),
                      sd   = sapply(daily[preds], sd))
for (p in preds) daily[[paste0(p, "_z")]] <- (daily[[p]] - scaling$mean[scaling$variable == p]) /
                                             scaling$sd[scaling$variable == p]

# ---- Step 7: station index (1, 2, 3) for the hierarchical model ------------
daily$station <- factor(daily$station)
daily$station_id <- as.integer(daily$station)

daily <- daily[order(daily$station, daily$date), ]

# ---- Step 8: quick sanity checks -------------------------------------------
cat("\nDays per station:\n"); print(table(daily$station))
cat("\nSummary of model variables:\n")
print(summary(daily[, c("log_pm25", "WSPM_z", "RAIN_z", "PRES_z")]))
cat("\nCorrelation among predictors (watch for very high values):\n")
print(round(cor(daily[, c("WSPM_z", "RAIN_z", "PRES_z")]), 2))

dir.create("figures", showWarnings = FALSE)
png("figures/eda_pm25.png", 1200, 450, res = 120)
par(mfrow = c(1, 2))
hist(daily$PM2.5, breaks = 60, main = "Daily PM2.5 (raw)", xlab = "ug/m^3", col = "grey80")
hist(daily$log_pm25, breaks = 60, main = "Daily log PM2.5", xlab = "log ug/m^3", col = "grey80")
dev.off()

# ---- Step 9: save ----------------------------------------------------------
write.csv(daily, file.path(data_dir, "air_daily.csv"), row.names = FALSE) # dataset the model will use 
write.csv(scaling, file.path(data_dir, "scaling.csv"), row.names = FALSE) # Since we standardize the predictors this is the csv to convert those untis back to their original ones
cat("\nSaved", nrow(daily), "rows to data/air_daily.csv\n")
