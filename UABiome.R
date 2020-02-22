#Microbiome Analysis Tool
#Shelby Nelson: Team 19060 Senior Design 
#2/9/2020

library(tidyverse)
library(tidyr)
library(dplyr)
library(readr)

#-----Reading in Files--------------------------------------------------------------

asv_table_f <- read_tsv("../Data/ASV_table.tsv")%>% rename("SampleID"="X1")
taxa_f <- read_tsv("../Data/taxa.tsv")
patient_f <- read_tsv("../Data/Patient_Data.txt")

#-----Combining Files into a Full Table---------------------------------------------
half_table <- inner_join(patient_f,                    #Joining Patient file and asv table by SampleID
                         asv_table_f, 
                         by = "SampleID") %>% 
  pivot_longer(c(ASV_1:ASV_421),                       #Pivot only the ASV numbers
               names_to = "ASV_num", 
               values_to = "count")

full_table <- inner_join(half_table,                   #Join the half table with the taxa file
                         taxa_f, 
                         by = "ASV_num")

full_table[is.na(full_table)] <- "NA"                  #Converting all blanks with NA
view(full_table)
#-----Save Full Text File----------------------------------------------------------
write_tsv(full_table,"C:\Users\jcmas\Documents\University of Arizona\Engineering\Spring 2020\ENGR 498A\Data\full_table.txt")

#-----Patient BT023 Samples--------------------------------------------------------
all_BT023_Samples <- filter(full_table,                #All Samples for this Patient
                            Patient == "BT023") 
single_BT023_Sample <- filter(full_table,              #Single sample for this patient with all ASvs= 0 removed
                              Patient == "BT023", 
                              SampleID == "Calton305", 
                              count != 0) 
view(single_BT023_Sample)
grouped_genus <- group_by(single_BT023_Sample, 
                          Genus, 
                          Day_of_Life) %>%
  mutate(Genus_Count = sum(count))                     #Sum each genus by the count

total_count <- sum(grouped_genus$count)                #sum the total count column

#genus_Frequency <- mutate(grouped_genus,               #creating a frequency column
#                          frequency = (Genus_Count/total_count)*100) %>% 
#  distinct(frequency)                                  #removing the duplicate percentages

#ggplot(data = genus_Frequency) +                       # column plot
#  geom_col(mapping = aes(x = Day_of_Life,
#                         y = frequency,
#                         fill = Genus))

uniq_genus<-grouped_genus%>%select(-ASV_num,-Sequence,-count) %>% unique()

grouped_genus %>% ggplot(aes(fill=Genus, y=Genus_Count, x=SampleID))+
  geom_bar(position = "fill", stat="identity")
