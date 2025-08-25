path <- 'C:/Users/EMCDUF01/Industrial Economics, Inc/IEc EPA CSIB Share - IEc EPA CSIB Share/03_Global_AQ_Study/Data/GCAP2 Raw Data 2025/CESM'

#base
path_full <- file.path(path, 'PM25.noDSTnoSSA.base.nc')
test <- nc_open(path_full)
pmbase <- ncvar_get(test,'PM25')

#126
path_full <- file.path(path, 'PM25.noDSTnoSSA.clim126.nc')
test <- nc_open(path_full)
pm126 <- ncvar_get(test,'PM25')

#245
path_full <- file.path(path, 'PM25.noDSTnoSSA.clim245.nc')
test <- nc_open(path_full)
pm245 <- ncvar_get(test,'PM25')

#585
path_full <- file.path(path, 'PM25.noDSTnoSSA.clim585.nc')
test <- nc_open(path_full)
pm585 <- ncvar_get(test,'PM25')

#calculate deltas
data<- rbind(data.frame(Values = as.table(pmbase-pmbase),Group = 'pmbase'),
             data.frame(Values = as.table(pm126-pmbase),Group = 'pm126'),
             data.frame(Values = as.table(pm245-pmbase), Group = 'pm245'),
             data.frame(Values = as.table(pm585-pmbase),Group = 'pm585'))

p<- ggplot(data, aes(Values.Freq, fill = Group)) + geom_histogram(position = "identity", alpha = 0.4, bins = 50)

print('PM2.5 Base (min, max, mean, median)')
print(min(data$Values.Freq[data$Group=='pmbase']))
print(max(data$Values.Freq[data$Group=='pmbase']))
print(mean(data$Values.Freq[data$Group=='pmbase']))
print(median(data$Values.Freq[data$Group=='pmbase']))

print('PM2.5 126 (min, max, mean, median)')
print(min(data$Values.Freq[data$Group=='pm126']))
print(max(data$Values.Freq[data$Group=='pm126']))
print(mean(data$Values.Freq[data$Group=='pm126']))
print(median(data$Values.Freq[data$Group=='pm126']))

print('PM2.5 245 (min, max, mean, median)')
print(min(data$Values.Freq[data$Group=='pm245']))
print(max(data$Values.Freq[data$Group=='pm245']))
print(mean(data$Values.Freq[data$Group=='pm245']))
print(median(data$Values.Freq[data$Group=='pm245']))

print('PM2.5 585 (min, max, mean, median)')
print(min(data$Values.Freq[data$Group=='pm585']))
print(max(data$Values.Freq[data$Group=='pm585']))
print(mean(data$Values.Freq[data$Group=='pm585']))
print(median(data$Values.Freq[data$Group=='pm585']))

#count number of rows with specific conditions
nrow(data[data$Group=='pm126' & data$Values.Freq < -10,])/259200

#p1<- ggplot(as.data.frame(as.table(pm126)), aes(x=Freq))+geom_histogram()
#p2 <- ggplot(as.data.frame(as.table(pm245)), aes(x=Freq))+geom_histogram()

#pbig <- p1+p2                        
