# Set working directory and load necessary packages
setwd("~/Box/Barott lab/Data/RTE2016/B-NB Temp")
library(readr)
library(lubridate)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)

# Import data
U22 <- read_csv("U22 Corrected.csv") 
Coconut <- read_csv("Mokuoloe1612480NOAA_July2015_2016.csv")
A12 <- read_csv("12A.csv")
B12 <- read_csv("12B.csv")

# Rename columns
names(U22) <- c("Date", "Time", "Time24", "InnerBay", "OuterBay")
names(Coconut) <- c("Date", "Time", "Time24", "InnerBay")
names(A12) <- c("Date", "Time", "Time24", "TempF", "TempC", "CorTempC")
names(B12) <- c("Date", "Time", "Time24", "TempF", "TempC", "CorTempC")

# Convert time to time format
U22$newTime <- as.POSIXct(U22$Time24, format = "%H:%M") 
U22$TimeFinal <- format(U22$newTime, "%H:%M") 
Coconut$newTime <- as.POSIXct(Coconut$Time24, format = "%H:%M") 
Coconut$TimeFinal <- format(Coconut$newTime, "%H:%M")
A12$newTime <- as.POSIXct(A12$Time24, format = "%H:%M") 
A12$TimeFinal <- format(A12$newTime, "%H:%M")
B12$newTime <- as.POSIXct(B12$Time24, format = "%H:%M") 
B12$TimeFinal <- format(B12$newTime, "%H:%M")

# Convert date to date format
U22$newDate <- as.Date(U22$Date, format = "%m/%d/%y")
Coconut$newDate <- as.Date(Coconut$Date, format = "%m/%d/%y")
A12$newDate <- as.Date(A12$Date, format = "%m/%d/%y")
B12$newDate <- as.Date(B12$Date, format = "%m/%d/%y")

# Merge date and time
U22$dttm <- as.POSIXct(paste(U22$newDate, U22$TimeFinal), format="%Y-%m-%d %H:%M") 
U22$dttm <- as_datetime(U22$dttm, tz = "America/New_York") 
Coconut$dttm <- as.POSIXct(paste(Coconut$newDate, Coconut$TimeFinal), format="%Y-%m-%d %H:%M") 
Coconut$dttm <- as_datetime(Coconut$dttm, tz = "America/New_York")  
A12$dttm <- as.POSIXct(paste(A12$newDate, A12$TimeFinal), format="%Y-%m-%d %H:%M") 
A12$dttm <- as_datetime(A12$dttm, tz = "America/New_York") 
B12$dttm <- as.POSIXct(paste(B12$newDate, B12$TimeFinal), format="%Y-%m-%d %H:%M") 
B12$dttm <- as_datetime(B12$dttm, tz = "America/New_York") 

# Assign group to continuous data
TimeDiff <- difftime(U22$dttm, lag(U22$dttm, default = U22$dttm[1]), units = "days")
U22$grp <- cumsum(ifelse(TimeDiff>2,1,0))
TimeDiff <- difftime(Coconut$dttm, lag(Coconut$dttm, default = Coconut$dttm[1]), units = "days")
Coconut$grp <- cumsum(ifelse(TimeDiff>2,1,0))
TimeDiff <- difftime(A12$dttm, lag(A12$dttm, default = A12$dttm[1]), units = "days")
A12$grp <- cumsum(ifelse(TimeDiff>2,1,0))
TimeDiff <- difftime(B12$dttm, lag(B12$dttm, default = B12$dttm[1]), units = "days")
B12$grp <- cumsum(ifelse(TimeDiff>2,1,0))

# Combine dataframes
AllLogs <- merge(merge(merge(Coconut, U22, by = "dttm", all = TRUE), A12, by = "dttm", all = TRUE), B12, by = "dttm", all = TRUE)
AllLogs <- AllLogs[ ,c(1,5,8:9,13:14,17:18,24,27:28,34,37:38)]
names(AllLogs) <- c("dttm", "PR1", "PR1Date", "PR1Group", "PR4", "PR13", "PR4Date", "PR4Group", "PR12", "PR12Date", "PR12Group", "B12", "B12Date", "B12Group")

# Structure data
Sub <- dplyr::select(AllLogs, dttm, PR1, PR4, PR12, PR13)
MeltSub <- melt(Sub, id.vars = c("dttm"))
MeltSub <- na.omit(MeltSub)

# Plot all hourly data
png("AllTemp.png", width = 200, height = 80, units = "mm", res = 500)
ggplot(MeltSub, aes(x = dttm, color = variable)) + 
  geom_line(aes(y = value), size = 0.35) +
  scale_color_manual(values = c("PR1" = "coral", "PR4" = "coral3", "PR12" = "skyblue2", "PR13" = "skyblue4")) +
  scale_x_datetime(labels = date_format("%b '%y"), date_breaks = "2 months") +
  theme(aspect.ratio = .3, axis.text=element_text(size=10, color = "black"), axis.title=element_text(size=12, color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.title = element_blank(), legend.key = element_blank(), 
        legend.text = element_text(size = 12, color = "black"), legend.position = "right", legend.key.size = unit(0.5, "cm")) +
  scale_fill_manual(values=alpha(c("coral", "coral3", "skyblue2", "skyblue4")), name = "Patch Reef") + 
  labs(y = "Temp (°C)", x = "Date", color = "Patch Reef", title = "")
dev.off()

# All data daily mean +/- max and min
U22Mean <- cbind(aggregate(AllLogs$PR4, by=AllLogs["PR4Date"], FUN=mean, na.rm = TRUE), aggregate(AllLogs$PR4, by=AllLogs["PR4Date"], FUN=max, na.rm = TRUE), aggregate(AllLogs$PR4, by=AllLogs["PR4Date"], FUN=min, na.rm = TRUE),
                 aggregate(AllLogs$PR13, by=AllLogs["PR4Date"], FUN=mean, na.rm = TRUE), aggregate(AllLogs$PR13, by=AllLogs["PR4Date"], FUN=max, na.rm = TRUE), aggregate(AllLogs$PR13, by=AllLogs["PR4Date"], FUN=min, na.rm = TRUE))
U22Mean <- U22Mean[ ,c(1:2,4,6,8,10,12)]
names(U22Mean) <- c("Date", "PR4", "PR4Max", "PR4Min", "PR13", "PR13Max", "PR13Min")
HIMBMean <- cbind(aggregate(AllLogs$PR1, by=AllLogs["PR1Date"], FUN=mean, na.rm = TRUE), aggregate(AllLogs$PR1, by=AllLogs["PR1Date"], FUN=max, na.rm = TRUE), aggregate(AllLogs$PR1, by=AllLogs["PR1Date"], FUN=min, na.rm = TRUE))
HIMBMean <- HIMBMean[ ,c(1:2,4,6)]
names(HIMBMean) <- c("Date", "PR1", "PR1Max", "PR1Min")
A12Mean <- cbind(aggregate(AllLogs$PR12, by=AllLogs["PR12Date"], FUN=mean, na.rm = TRUE), aggregate(AllLogs$PR12, by=AllLogs["PR12Date"], FUN=max, na.rm = TRUE), aggregate(AllLogs$PR12, by=AllLogs["PR12Date"], FUN=min, na.rm = TRUE))
A12Mean <- A12Mean[ ,c(1:2,4,6)]
names(A12Mean) <- c("Date", "PR12", "PR12Max", "PR12Min")
AllMean <- merge(merge(HIMBMean, U22Mean, by = "Date", all = TRUE), A12Mean, by = "Date", all = TRUE)
is.na(AllMean) <- do.call(cbind,lapply(AllMean, is.infinite))

# Structure data
MeanSub <- dplyr::select(AllMean, Date, PR1, PR4, PR12, PR13)
MeanMelt <- melt(MeanSub, id.vars = c("Date"))
MeanMelt <- na.omit(MeanMelt)

MaxSub <- dplyr::select(AllMean, Date, PR1Max, PR4Max, PR12Max, PR13Max)
MaxMelt <- melt(MaxSub, id.vars = c("Date"))
MaxMelt <- na.omit(MaxMelt)

MinSub <- dplyr::select(AllMean, Date, PR1Min, PR4Min, PR12Min, PR13Min)
MinMelt <- melt(MinSub, id.vars = c("Date"))
MinMelt <- na.omit(MinMelt)

MeanMelt <- cbind(MeanMelt, MaxMelt, MinMelt)
MeanMelt <- MeanMelt[ ,c(1:3,6,9)]
names(MeanMelt) <- c("Date", "PR", "Mean", "Max", "Min")

# Ignore PR1/PR12 after 10-1 
MeanMelt <- subset(MeanMelt, Date >= "2015-07-01")
MeanMelt$Mean <- ifelse(MeanMelt$Date > "2015-10-01" & MeanMelt$PR == "PR1", NA, MeanMelt$Mean)
MeanMelt$Mean <- ifelse(MeanMelt$Date > "2015-10-01" & MeanMelt$PR == "PR12", NA, MeanMelt$Mean)
MeanMelt$Max <- ifelse(MeanMelt$Date > "2015-10-01" & MeanMelt$PR == "PR1", NA, MeanMelt$Max)
MeanMelt$Max <- ifelse(MeanMelt$Date > "2015-10-01" & MeanMelt$PR == "PR12", NA, MeanMelt$Max)
MeanMelt$Min <- ifelse(MeanMelt$Date > "2015-10-01" & MeanMelt$PR == "PR1", NA, MeanMelt$Min)
MeanMelt$Min <- ifelse(MeanMelt$Date > "2015-10-01" & MeanMelt$PR == "PR12", NA, MeanMelt$Min)

# Plot Daily Mean
png("MeanTemp.png", width = 200, height = 80, units = "mm", res = 500)
ggplot(MeanMelt, aes(x = Date, color = PR)) + 
  geom_hline(yintercept = 28, linetype = "dashed") +
  geom_line(aes(y = Mean), size = 0.5) +
  geom_ribbon(aes(ymin = Min, ymax = Max, fill = PR), alpha = 0.2, linetype = "blank") +
  scale_color_manual(values = c("PR1" = "coral3", "PR4" = "coral", "PR12" = "skyblue4", "PR13" = "skyblue2")) +
  scale_x_date(labels = date_format("%b '%y"), date_breaks = "2 months", limits = as.Date(c("2015-06-20", "2016-10-06"))) + 
  theme(aspect.ratio = .3, axis.title=element_text(size=12, color = "black", face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.title = element_text(size = 10, color = "black", face = "bold"), legend.key = element_blank(), 
        legend.text = element_text(size = 10, color = "black"), legend.position = c(0.065, 0.2), legend.key.size = unit(0.3, "cm")) +
  scale_fill_manual(values=alpha(c("coral3", "coral", "skyblue4", "skyblue2")), name = "Patch Reef") + 
  labs(y = "Temperature (°C)", x = "Date", color = "Patch Reef", title = "") +
  guides(color=guide_legend(override.aes=list(fill=NA)))
dev.off()

# plot Daily Max
png("MaxTemp.png", width = 200, height = 80, units = "mm", res = 500)
ggplot(MeanMelt, aes(x = Date, color = PR)) + 
  geom_hline(yintercept = 28, linetype = "dashed") +
  geom_hline(yintercept = 29) +
  geom_line(aes(y = Max), size = 0.5) +
  scale_color_manual(values = c("PR1" = "coral3", "PR4" = "coral", "PR12" = "skyblue4", "PR13" = "skyblue2")) +
  scale_x_date(labels = date_format("%b '%y"), date_breaks = "2 months", limits = as.Date(c("2015-06-20", "2016-10-06"))) + 
  scale_y_continuous(breaks = c(24,26,28,30)) +
  theme(aspect.ratio = .3, axis.title=element_text(size=12, color = "black", face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.title = element_text(size = 10, color = "black", face = "bold"), legend.key = element_blank(), 
        legend.text = element_text(size = 10, color = "black"), legend.position = c(0.065, 0.2), legend.key.size = unit(0.3, "cm")) +
  scale_fill_manual(values=alpha(c("coral3", "coral", "skyblue4", "skyblue2")), name = "Patch Reef") + 
  labs(y = "Temperature (°C)", x = "Date", color = "Patch Reef", title = "") +
  guides(color=guide_legend(override.aes=list(fill=NA)))
dev.off()

##### Stats ##### 
## BLEACH 
BleachMean <- subset(MeanMelt, Date <= "2015-10-01" & Date >= "2015-07-01")
BleachMean$Range <- BleachMean$Max - BleachMean$Min

BleachStats <- cbind(aggregate(BleachMean$Mean, by = BleachMean["PR"], FUN = mean, na.rm = TRUE),
                     aggregate(BleachMean$Max, by = BleachMean["PR"], FUN = mean, na.rm = TRUE),
                     aggregate(BleachMean$Min, by = BleachMean["PR"], FUN = mean, na.rm = TRUE),
                     aggregate(BleachMean$Mean, by = BleachMean["PR"], FUN = sd, na.rm = TRUE),
                     aggregate(BleachMean$Range, by = BleachMean["PR"], FUN = mean, na.rm = TRUE))
BleachStats <- BleachStats[, c(1:2,4,6,8,10)]
names(BleachStats) <- c("Reef", "Mean", "Max", "Min", "SD", "Range")

sum(BleachMean$Max >= 29 & BleachMean$PR == "PR1") # 46
sum(BleachMean$Max >= 29 & BleachMean$PR == "PR12") # 24

## RECOVERY 
Recov <- subset(MeanMelt, Date >= "2015-11-01" & Date <= "2016-02-01")
Recov$Range <- Recov$Max - Recov$Min

RecovStats <- cbind(aggregate(Recov$Mean, by = Recov["PR"], FUN = mean, na.rm = TRUE),
                     aggregate(Recov$Max, by = Recov["PR"], FUN = mean, na.rm = TRUE),
                     aggregate(Recov$Min, by = Recov["PR"], FUN = mean, na.rm = TRUE),
                     aggregate(Recov$Mean, by = Recov["PR"], FUN = sd, na.rm = TRUE),
                     aggregate(Recov$Range, by = Recov["PR"], FUN = mean, na.rm = TRUE))
RecovStats <- RecovStats[, c(1:2,4,6,8,10)]
names(RecovStats) <- c("Reef", "Mean", "Max", "Min", "SD", "Range")

# DHW
## PR1
Bleach <- subset(AllLogs, dttm <= "2015-10-01" & dttm >= "2015-07-01" & PR1 != "NA")
Bleach$Stress <- Bleach$PR1-28
DHWPR128 <- subset(Bleach, Stress >= 1)
DHWPR128$Stress <- floor(DHWPR128$Stress)
sum(DHWPR128$Stress)/168 # 5.839286 

## PR12
Bleach$Stress3 <- Bleach$PR12-28
DHWPR1228 <- subset(Bleach, Stress3 >= 1)
DHWPR1228$Stress3 <- floor(DHWPR1228$Stress3)
sum(DHWPR1228$Stress3)/168 # 1.77381 

# Mean daily max - Bleach
Bleach <- subset(AllMean, Date <= "2015-10-01" & Date >= "2015-07-01" & PR1 != "NA")
shapiro.test(Bleach$PR1Max) # p = 0.00087 < 0.05 -> non-normal
shapiro.test(Bleach$PR12Max) # p = 0.003879 < 0.05 -> non-normal
var.test(Bleach$PR1Max, Bleach$PR12Max) # p = 0.9523 > 0.05 -> equal variances
t.test(Bleach$PR1Max, Bleach$PR12Max, var.equal = T) # p = 4.599e-05 -> sig 
wilcox.test(Bleach$PR1Max, Bleach$PR12Max) # p = 1.477e-05 -> sig

# Mean daily min - Bleach
shapiro.test(Bleach$PR1Min) # p = 0.006793 < 0.05 -> non-normal
shapiro.test(Bleach$PR12Min) # p = 5.222e-05 < 0.05 -> non-normal
var.test(Bleach$PR1Min, Bleach$PR12Min) # p = .1464 > 0.05 -> equal variances
t.test(Bleach$PR1Min, Bleach$PR12Min, var.equal = T) # p = 3.828e-07 -> sig 
wilcox.test(Bleach$PR1Min, Bleach$PR12Min) # p = 5.449e-08 -> sig

# Mean daily range - Bleach
Bleach$PR1Range <- Bleach$PR1Max - Bleach$PR1Min
shapiro.test(Bleach$PR1Range) # p = 0.01125 > 0.05 -> normal
Bleach$PR12Range <- Bleach$PR12Max - Bleach$PR12Min
shapiro.test(Bleach$PR12Range) # p = 0.1532 > 0.05 -> normal
var.test(Bleach$PR1Range, Bleach$PR12Range) # p = 0.001481 < 0.05 -> non-equal variances
t.test(Bleach$PR1Range, Bleach$PR12Range, var.equal = F) # p = 0.03722 -> not sig 
wilcox.test(Bleach$PR1Range, Bleach$PR12Range) # p = 0.1188 -> not sig

# Mean daily max - Recov
Recov <- subset(AllMean, Date >= "2015-11-01" & Date <= "2016-02-01")
shapiro.test(Recov$PR1Max) # p = 5.508e-05 < 0.05 -> non-normal
shapiro.test(Recov$PR12Max) # p = 2.77e-05 < 0.05 -> non-normal
var.test(Recov$PR1Max, Recov$PR12Max) # p = 0.162 > 0.05 -> equal variances
t.test(Recov$PR1Max, Recov$PR12Max, var.equal = T) # p = 0.04698 -> not sig 
wilcox.test(Recov$PR1Max, Recov$PR12Max) # p = 0.02185 -> not sig
shapiro.test(Recov$PR4Max) # p = 0.001465 < 0.05 -> non-normal
shapiro.test(Recov$PR13Max) # p = 0.0008875 < 0.05 -> non-normal
var.test(Recov$PR4Max, Recov$PR13Max) # p = 0.02918 > 0.05 -> equal variances
t.test(Recov$PR4Max, Recov$PR13Max, var.equal = T) # p = 0.8119 -> not sig 
wilcox.test(Recov$PR4Max, Recov$PR13Max) # p = 0.5906 -> not sig

# Mean daily min - Recov
shapiro.test(Recov$PR1Min) # p = 0.000138 < 0.05 -> non-normal
shapiro.test(Recov$PR12Min) # p = 5.078e-05 < 0.05 -> non-normal
var.test(Recov$PR1Min, Recov$PR12Min) # p = 0.3962 > 0.05 -> equal variances
t.test(Recov$PR1Min, Recov$PR12Min, var.equal = T) # p = 0.04819 -> not sig 
wilcox.test(Recov$PR1Min, Recov$PR12Min) # p = 0.01969 -> not sig
shapiro.test(Recov$PR4Min) # p = 0.0005599 < 0.05 -> non-normal
shapiro.test(Recov$PR13Min) # p = 0.000358 < 0.05 -> non-normal
var.test(Recov$PR4Min, Recov$PR13Min) # p = 0.1416 > 0.05 -> equal variances
t.test(Recov$PR4Min, Recov$PR13Min, var.equal = T) # p = 0.2608 -> not sig 
wilcox.test(Recov$PR4Min, Recov$PR13Min) # p = 0.2152 -> not sig

# Mean daily range - Recov
Recov$PR1Range <- Recov$PR1Max - Recov$PR1Min
shapiro.test(Recov$PR1Range) # p = 0.09899 > 0.05 -> normal
Recov$PR12Range <- Recov$PR12Max - Recov$PR12Min
shapiro.test(Recov$PR12Range) # p = 0.0005528 < 0.05 -> non-normal
var.test(Recov$PR1Range, Recov$PR12Range) # p = 0.03293 > 0.05 -> equal variances
t.test(Recov$PR1Range, Recov$PR12Range, var.equal = T) # p = 0.9935 -> not sig 
wilcox.test(Recov$PR1Range, Recov$PR12Range) # p = 0.4455 -> not sig
Recov$PR4Range <- Recov$PR4Max - Recov$PR4Min
shapiro.test(Recov$PR4Range) # p = 0.0007276 < 0.05 -> non-normal
Recov$PR13Range <- Recov$PR13Max - Recov$PR13Min
shapiro.test(Recov$PR13Range) # p = 0.04107 < 0.05 -> non-normal
var.test(Recov$PR4Range, Recov$PR13Range) # p = 0.0001996 < 0.05 -> non-equal variances
t.test(Recov$PR4Range, Recov$PR13Range, var.equal = F) # p = 0.0.0103 -> not sig 
wilcox.test(Recov$PR4Range, Recov$PR13Range) # p = 0.01734 -> not sig

# Extract overlapping dates
png("scattercompare.png", width = 200, height = 200, units = "mm", res = 500)
par(mfrow = c(2,2))
OverlapInner <- subset(AllMean, (!is.na(AllMean$PR1Max)) & (!is.na(AllMean$PR4Max)))
plot(OverlapInner$PR1Max, OverlapInner$PR4Max, ylim = c(23, 29), xlim = c(23, 29), ylab = "PR4 Daily Max Temp (°C)", xlab = "PR1 Daily Max Temp (°C)")
abline(a = 0, b = 1)
fit <- lm(PR1Max ~ PR4Max, data = OverlapInner)
legend(22.5, 28.3, bty = "n", legend = paste("R2:",format(summary(fit)$adj.r.squared, digits=4)))
cf <- round(coef(fit), 2) 
eq <- paste0("y = ", cf[1],
             ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x")
mtext(eq, 3, line = -2, adj = 0.08)
OverlapOuter <- subset(AllMean, (!is.na(AllMean$PR12Max)) & (!is.na(AllMean$PR13Max)))
plot(OverlapInner$PR12Max, OverlapInner$PR13Max, ylim = c(23, 29), xlim = c(23, 29), ylab = "PR13 Daily Max Temp (°C)", xlab = "PR12 Daily Max Temp (°C)")
abline(a = 0, b = 1)
fit <- lm(PR12Max ~ PR13Max, data = OverOuter)
legend(22.5, 28.3, bty = "n", legend = paste("R2:",format(summary(fit)$adj.r.squared, digits=4)))
cf <- round(coef(fit), 2) 
eq <- paste0("y = ", cf[1],
             ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x")
mtext(eq, 3, line = -2, adj = 0.08)
plot(OverlapInner$PR1, OverlapInner$PR4, ylim = c(23, 29), xlim = c(23, 29), ylab = "PR4 Daily Mean Temp (°C)", xlab = "PR1 Daily Mean Temp (°C)")
abline(a = 0, b = 1)
fit <- lm(PR1 ~ PR4, data = OverlapInner)
legend(22.5, 28.3, bty = "n", legend = paste("R2:",format(summary(fit)$adj.r.squared, digits=4)))
cf <- round(coef(fit), 2) 
eq <- paste0("y = ", cf[1],
             ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x")
mtext(eq, 3, line = -2, adj = 0.08)
plot(OverlapInner$PR12, OverlapInner$PR13, ylim = c(23, 29), xlim = c(23, 29), ylab = "PR13 Daily Mean Temp (°C)", xlab = "PR12 Daily Mean Temp (°C)")
abline(a = 0, b = 1)
fit <- lm(PR12 ~ PR13, data = OverOuter)
legend(22.5, 28.3, bty = "n", legend = paste("R2:",format(summary(fit)$adj.r.squared, digits=4)))
cf <- round(coef(fit), 2) 
eq <- paste0("y = ", cf[1],
             ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x")
mtext(eq, 3, line = -2, adj = 0.08)
dev.off()



