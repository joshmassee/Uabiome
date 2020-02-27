#Microbiome Analysis Tool
#Team 19060 Senior Design 
#2/9/2020

library(tidyverse)
library(tidyr)
library(dplyr)
library(vegan)


#-----Reading in Files--------------------------------------------------------------
asv_table_f <- read_tsv("C:/Users/Shelbeezy/work/uabiome_data/ASV_table.tsv")
taxa_f <- read_tsv("C:/Users/Shelbeezy/work/uabiome_data/taxa.tsv")
patient_f <- read_tsv("C:/Users/Shelbeezy/work/uabiome_data/Patient_Data.txt")

#-----Combining Files into a Full Table---------------------------------------------
half_table <- inner_join(patient_f,                    #Joining Patient file and asv table by SampleID
                         asv_table_f, 
                         by = "SampleID") %>% 
  pivot_longer(c(ASV_1:ASV_421),                       #Pivot only the ASV numbers
               names_to = "ASV_num", 
               values_to = "count")

full_table <- inner_join(half_table,                   #Join the half table with the taxa file
                         taxa_f, 
                         by = "ASV_num") %>% 
  filter(count != 0,                                   #Filter out counts that are equal to zero and Kingdoms that are not Bacteria
         Kingdom == "Bacteria") %>%
  group_by(SampleID, 
           Genus) %>%
  mutate(Genus_Count = sum(count)) %>% 
  select(-ASV_num,
         -Sequence,
         -count) %>%
  unique()

full_table[is.na(full_table)] <- "Other"                  #Converting all blanks with Other

view(full_table)
#-----------Table of genus averages across all samples and patients----------------------
average_sum_full <- full_table %>%                          #Finding sum of genus count
  group_by(Genus) %>% 
  mutate(Genus_Count = sum(Genus_Count)) %>%
  select(Genus, 
         Genus_Count) %>%
  unique()

average_tally_full <- full_table %>%                        #Finding total number of that genus
  group_by(Genus) %>% 
  tally() 

average_genus_full <- inner_join(average_sum_full,              #Finding Average of each Genus type
                            average_tally_full, 
                            by = "Genus") %>%
  mutate(Average_Count = Genus_Count/n) %>% 
  select(Genus, 
         Average_Count)

view(average_genus_full)

#-----------Single Patient Information-----------------------------------------------
sample_1 <- filter(full_table,                   
                   Patient == "BT023",
                   SampleID == "Calton301") %>%
  mutate(Genus_Ratio = Genus_Count/sum(sample_1$Genus_Count))
view(sample_1)
sample_2 <- filter(full_table,                   
                   Patient == "BT023",
                   SampleID == "Calton302")
sample_3 <- filter(full_table,                   
                   Patient == "BT023",
                   SampleID == "Calton303")
sample_4 <- filter(full_table,                   
                   Patient == "BT023",
                   SampleID == "Calton304")
sample_5 <- filter(full_table,                   
                   Patient == "BT023",
                   SampleID == "Calton305")

desc_sample_1 <- sample_1[order(desc(sample_1$Genus_Count)),]
desc_sample_2 <- sample_2[order(desc(sample_2$Genus_Count)),]
desc_sample_3 <- sample_3[order(desc(sample_3$Genus_Count)),]
desc_sample_4 <- sample_4[order(desc(sample_4$Genus_Count)),]
desc_sample_5 <- sample_5[order(desc(sample_5$Genus_Count)),]

patient_1 <- rbind(desc_sample_1,                      #All of patient 1 in descending order
                   desc_sample_2,
                   desc_sample_3,
                   desc_sample_4,
                   desc_sample_5)

view(patient_1)

top10_each <- patient_1 %>%                           # Top 10 of each sample of patient 1 in one dataframe for plotting
  group_by(SampleID) %>%
  top_n(n = 10, 
        Genus_Count)

view(top10_each)

write_tsv(top10_each, "C:/Users/Shelbeezy/work/uabiome_data/top10_PatientBT023.txt")

#-------------Averages in progress------------------------------------------------------
# average_sum <- patient_1 %>%                          #Finding sum of genus count
#   group_by(Genus) %>% 
#   mutate(Genus_Count = sum(Genus_Count)) %>%
#   select(Genus, 
#          Genus_Count) %>%
#   unique()
# 
# average_tally <- patient_1 %>%                        #Finding total number of that genus
#   group_by(Genus) %>% 
#   tally() 
# 
# average_genus <- inner_join(average_sum,              #Finding Average of each Genus type
#                             average_tally, 
#                             by = "Genus") %>%
#   mutate(Average_Count = Genus_Count/n) %>% 
#   select(Genus, 
#          Average_Count)
# 
# view(average_genus)

#---------Alpha Diversity in Progress----------------------------------------------------

#alpha <- diversity(sample_1$Genus_Ratio)
#view(alpha)
#view(sample_1)
#---------Single Plotting---------------------------------------------------------------

#---------Single: Plotting Top 10 Most Abundant for Sample 1-----------------------------
top10_each %>% 
  filter(SampleID == "Calton301") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID)) +                             
  geom_bar(position = "fill", 
           stat = "identity") +
  labs(title = "Total Genus Abundance for Sample 1",
       x = "Sample ID",
       y = "Frequency of Abundance",
       fill = "Genus")

#---------Single: Plotting Top 10 Most Abundant for Sample 2-----------------------------
top10_each %>% 
  filter(SampleID == "Calton302") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID)) +                             
  geom_bar(position = "fill", 
           stat = "identity") +
  labs(title = "Total Genus Abundance for Sample 2",
       x = "Sample ID",
       y = "Frequency of Abundance",
       fill = "Genus")

#---------Single: Plotting Top 10 Most Abundant for Sample 3-----------------------------
top10_each %>% 
  filter(SampleID == "Calton303") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID)) +                             
  geom_bar(position = "fill", 
           stat = "identity") +
  labs(title = "Total Genus Abundance for Sample 3",
       x = "Sample ID",
       y = "Frequency of Abundance",
       fill = "Genus")

#---------Single: Plotting Top 10 Most Abundant for Sample 4-----------------------------
top10_each %>% 
  filter(SampleID == "Calton304") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID)) +                             
  geom_bar(position = "fill", 
           stat = "identity") +
  labs(title = "Total Genus Abundance for Sample 4",
       x = "Sample ID",
       y = "Frequency of Abundance",
       fill = "Genus")

#---------Single: Plotting Top 10 Most Abundant for Sample 5-----------------------------
top10_each %>% 
  filter(SampleID == "Calton305") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID)) +                             
  geom_bar(position = "fill", 
           stat = "identity") +
  labs(title = "Total Genus Abundance for Sample 5",
       x = "Sample ID",
       y = "Frequency of Abundance",
       fill = "Genus")

#---------Temporal Plotting---------------------------------------------------------------

#---------Temporal: Plotting Top 10 Most Abundant by SampleID-----------------------------
top10_each %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID)) +                             #This shows a condensed temporal view
  geom_bar(position = "fill", 
           stat = "identity") +
  labs(title = "Total Genus Abundance per Sample",
       subtitle = "Condensed",
       x = "Sample ID",
       y = "Frequency of Abundance",
       fill = "Genus")

#--------Temporal: Plotting Top 10 Most Abundant by Day of Life---------------------------
top10_each %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = Day_of_Life)) +                          #This shows a temporal view seperated by time
  geom_bar(position = "fill", 
           stat = "identity") +
  labs(title = "Total Genus Abundance per Sample",
       subtitle = "Spaced out over Time",
       x = "Date",
       y = "Frequency of Abundance",
       fill = "Genus")

