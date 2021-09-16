#Install and load dismo
install.packages("dismo")
library(dismo)
library(dplyr)
library(ggplot2)
#Extract all data for blacknose dace
LMB <- gbif(genus = "Micropterus", species = "salmoides")
CIS <- gbif(genus = "Coregonus", species = "artedi")
WAE <- gbif(genus = "Sander", species = "vitreus")

#### LMB DATA #### 
#Extract sampling columns "sampleSizeUnit", "sampleSizeValue", "samplingProtocol", "samplingEffort"
#depending on species these have different combos of these columns
LMB_dat <- select(LMB, adm1, adm2, basisOfRecord, cloc, datasetName, hostingOrganizationKey, institutionCode, institutionID, institutionKey, ownerInstitutionCode,  sampleSizeUnit, sampleSizeValue ,  samplingProtocol, samplingEffort, lat, lon, year )
LMB_dat$basisOfRecord<-as.factor(LMB_dat$basisOfRecord) #99 fossil and 12617 preserved from 99,900 records
ggplot(LMB_dat, aes(x=basisOfRecord)) +
  geom_bar() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45))

#look at a table of where the record came from 
sum_table<-LMB_dat %>%
  count(basisOfRecord) 

#Remove rows where there's no sampling info
sampling <- select(LMB_dat, basisOfRecord,sampleSizeUnit, sampleSizeValue ,  samplingProtocol, samplingEffort, year )
sampling$basisOfRecord<-as.factor(sampling$basisOfRecord) 
complete_sampling<-sampling[complete.cases(sampling[ , 2:5]),]
summary(complete_sampling$basisOfRecord)

#### CIS  DATA #### 
CIS_dat <- select(CIS, adm1, adm2, basisOfRecord, cloc, datasetName, hostingOrganizationKey, institutionCode, institutionID, institutionKey, ownerInstitutionCode, sampleSizeValue ,  samplingProtocol, samplingEffort, lat, lon, year )
CIS_dat$basisOfRecord<-as.factor(CIS_dat$basisOfRecord) #99 fossil and 12617 preserved from 99,900 records
ggplot(CIS_dat, aes(x=basisOfRecord)) +
  geom_bar() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45))

#look at a table of where the record came from 
sum_table<-CIS_dat %>%
  count(basisOfRecord) 

#Remove rows where there's no sampling info
sampling <- select(CIS_dat, basisOfRecord, sampleSizeValue ,  samplingProtocol, samplingEffort, year )
sampling$basisOfRecord<-as.factor(sampling$basisOfRecord) 
sampling_protocal<-sampling[complete.cases(sampling[ , 3]),]
sampling_effort<-sampling[complete.cases(sampling[ , 4]),]
summary(sampling_protocal$basisOfRecord)
summary(sampling_effort$basisOfRecord)

#### CIS  DATA #### 
WAE_dat <- select(WAE, adm1, adm2, basisOfRecord, cloc, datasetName, hostingOrganizationKey, institutionCode, institutionID, institutionKey, ownerInstitutionCode, samplingProtocol, samplingEffort, lat, lon, year )
WAE_dat$basisOfRecord<-as.factor(WAE_dat$basisOfRecord) #99 fossil and 12617 preserved from 99,900 records
ggplot(WAE_dat, aes(x=basisOfRecord)) +
  geom_bar() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45))
#look at a table of where the record came from 
sum_table<-WAE_dat %>%
  count(basisOfRecord) 

#Remove rows where there's no sampling info
sampling <- select(WAE_dat, basisOfRecord,  samplingProtocol, samplingEffort, year )
sampling$basisOfRecord<-as.factor(sampling$basisOfRecord) 
sampling_protocal<-sampling[complete.cases(sampling[ , 2]),]
sampling_effort<-sampling[complete.cases(sampling[ , 3]),]
summary(sampling_protocal$basisOfRecord)
summary(sampling_effort$basisOfRecord)

#### UMMZ data from iDigBio #### 

ummz<-read.csv('/Users/katelynking/Desktop/ummz_dat/occurrence_raw.csv')
ummz_sampling<-select(ummz, dwc.basisOfRecord, dwc.continent, dwc.country, dwc.family, dwc.genus, dwc.specificEpithet, dwc.sampleSizeValue, dwc.samplingEffort, dwc.samplingProtocol, dwc.fieldNotes, dwc.habitat, dwc.lifeStage, dwc.occurrenceDetails, obis.ExtendedMeasurementOrFact, dwc.year, dwc.waterBody, dwc.decimalLongitude, dwc.decimalLatitude)
sampling_protocal<-filter(ummz_sampling, dwc.samplingProtocol != "") # overall 72% have a gear

LMB_ummz <- filter(ummz_sampling,   dwc.genus == "Micropterus" & dwc.specificEpithet == "salmoides") #1759 records
sampling_protocal<-filter(LMB_ummz, dwc.samplingProtocol != "")

CIS_ummz <- filter(ummz_sampling,   dwc.genus == "Coregonus" & dwc.specificEpithet == "artedi") #1190
sampling_protocal<-filter(CIS_ummz, dwc.samplingProtocol != "")

WAE_ummz <- filter(ummz_sampling,   dwc.genus == "Sander" & dwc.specificEpithet == "vitreus") #345
sampling_protocal<-filter(WAE_ummz, dwc.samplingProtocol != "")
