#Microbiome Analysis Tool
#Team 19060 Senior Design 
#2/9/2020

library(tidyverse)
library(tidyr)
library(dplyr)
library(vegan)
library(ggpubr)
library(scales)
library(treemap)
library(waffle)

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

proportion_table <- full_table %>%                          
  group_by(SampleID) %>% 
  mutate(total_count = sum(Genus_Count)) %>% 
  ungroup() %>%
  unique()  %>% 
  mutate(proportion = Genus_Count/total_count)

view(proportion_table)

#-----------Single Patient Information-----------------------------------------------
ordered <- proportion_table %>% 
  filter(Patient == "A") %>%
  group_by(SampleID) %>% 
  arrange(desc(Genus_Count), 
          .by_group = TRUE)

view(ordered)

write_tsv(ordered, "C:/Users/Shelbeezy/work/uabiome_data/Patient_A.txt")

top10_each <- ordered %>%                           # Top 10 of each sample of patient 1 in one dataframe for plotting
  group_by(SampleID) %>%
  top_n(n = 10, 
        Genus_Count)

view(top10_each)

write_tsv(top10_each, "C:/Users/Shelbeezy/work/uabiome_data/top10_Patient_A.txt")

#---------Single: Plotting Top 10 Most Abundant for Sample 1-----------------------------
top10_each %>% 
  filter(SampleID == "Sample 1") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID),
         width = 1) +                             
  geom_bar(position = "fill", 
           stat = "identity",
           width = 0.6) + 
  scale_y_continuous(labels = percent) +
  labs(title = "Total Genus Abundance for Sample 1",
       x = "Sample Number",
       y = "Frequency of Abundance",
       fill = "Genus") + 
  theme(plot.title = element_text(size=18, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.text = element_text(size =11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face="bold"))

#---------Single: Plotting Top 10 Most Abundant for Sample 2-----------------------------
top10_each %>% 
  filter(SampleID == "Sample 2") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID),
         width = 1) +                             
  geom_bar(position = "fill", 
           stat = "identity",
           width = 0.6) + 
  scale_y_continuous(labels = percent) +
  labs(title = "Total Genus Abundance for Sample 2",
       x = "Sample Number",
       y = "Frequency of Abundance",
       fill = "Genus") + 
  theme(plot.title = element_text(size=18, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.text = element_text(size =11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face="bold"))

#---------Single: Plotting Top 10 Most Abundant for Sample 3-----------------------------
top10_each %>% 
  filter(SampleID == "Sample 3") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID),
         width = 1) +                             
  geom_bar(position = "fill", 
           stat = "identity",
           width = 0.6) + 
  scale_y_continuous(labels = percent) +
  labs(title = "Total Genus Abundance for Sample 3",
       x = "Sample Number",
       y = "Frequency of Abundance",
       fill = "Genus") + 
  theme(plot.title = element_text(size=18, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.text = element_text(size =11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face="bold"))

#---------Single: Plotting Top 10 Most Abundant for Sample 4-----------------------------
top10_each %>% 
  filter(SampleID == "Sample 4") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID),
         width = 1) +                             
  geom_bar(position = "fill", 
           stat = "identity",
           width = 0.6) + 
  scale_y_continuous(labels = percent) +
  labs(title = "Total Genus Abundance for Sample 4",
       x = "Sample Number",
       y = "Frequency of Abundance",
       fill = "Genus") + 
  theme(plot.title = element_text(size=18, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.text = element_text(size =11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face="bold"))

#---------Single: Plotting Top 10 Most Abundant for Sample 5-----------------------------
top10_each %>% 
  filter(SampleID == "Sample 5") %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID),
         width = 1) +                             
  geom_bar(position = "fill", 
           stat = "identity",
           width = 0.6) + 
  scale_y_continuous(labels = percent) +
  labs(title = "Total Genus Abundance for Sample 5",
       x = "Sample Number",
       y = "Frequency of Abundance",
       fill = "Genus") + 
  theme(plot.title = element_text(size=18, face="bold.italic"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.text = element_text(size =11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face="bold"))

#--------------Testing Plots--------------------------------------------

top10_each %>%                                   #Pie Chart
  filter(SampleID == "Sample 2") %>%
  ggplot(aes(x="", y=proportion, fill=Genus)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_void()+
  labs(title = "Total Genus Abundance for Sample 2",
       fill = "Genus") +
  theme(plot.title = element_text(size=18, face="bold.italic", hjust = .65, vjust = -5),
        legend.text = element_text(size =11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"))

top10_each %>%                                 #Treemap
  filter(SampleID == "Sample 4") %>%
  treemap(index = "Genus",
          vSize="Genus_Count",
          type="index",
          title = "Total Genus Abundance in Sample 4",
          fontsize.title = 18)

top10_each %>%
  ggplot(aes(x = Genus, 
             y = proportion * 100, 
             color = SampleID))+
  geom_point(size = 5, 
             alpha = 0.6) + 
  labs(title = "Total Genus Abundance Across All Samples",
       x = "Genus",
       y = "Frequency of Abundance") +
  scale_color_discrete(name = "Sample Number") +
  theme(axis.text.x = element_text(angle = 55, 
                                   hjust = 1),
        plot.title = element_text(size = 18, face = "bold.italic"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face = "bold"))

top10_each %>%
  ggplot(aes(x = Genus, 
             y = proportion * 100, 
             color = SampleID))+
  geom_point(size = 5, 
             alpha = 0.6) + 
  labs(title = "Total Genus Abundance Across All Samples",
       x = "Genus",
       y = "Frequency of Abundance") +
  scale_color_discrete(name = "Sample Number") +
  theme(axis.text.x = element_text(angle = 55, 
                                 hjust = 1),
        plot.title = element_text(size = 18, face = "bold.italic"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face = "bold")) +
  coord_flip()
 
top10_each %>%
  ggplot(aes(x = SampleID, 
             y = Genus)) +
  geom_point(aes(color = proportion * 100), 
             size = 10) + 
  labs(title = "Total Genus Abundance Across All Samples",
       x = "Sample Number",
       y = "Genus",
       fill = "Genus") +
  scale_colour_gradient(low = "#1ee6fc", 
                        high = "Black", 
                        name = "Frequency of\n Abundance") +
  theme(plot.title = element_text(size = 18, face = "bold.italic"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face = "bold"))

top10_each %>%
  ggplot(aes(x = SampleID, 
             y = Genus)) +
  geom_point(aes(color = proportion * 100), 
             size = 10) + 
  labs(title = "Total Genus Abundance Across All Samples",
       x = "Sample Number",
       y = "Genus",
       fill = "Genus") +
  scale_colour_gradient(low = "orange", 
                        high = "purple", 
                        name = "Frequency of\n Abundance") +
  theme(plot.title = element_text(size = 18, face = "bold.italic"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face = "bold"))


#---------Temporal Plotting---------------------------------------------------------------

#---------Temporal: Plotting Top 10 Most Abundant by SampleID-----------------------------
top10_each %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = SampleID)) +                             #This shows a condensed temporal view
  geom_bar(position = "fill", 
           stat = "identity") +
  labs(title = "Total Genus Abundance per Sample",
       x = "Sample Number",
       y = "Frequency of Abundance",
       fill = "Genus") +
  scale_y_continuous(labels = percent) +
  theme(plot.title = element_text(size = 18, face = "bold.italic"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11, face = "bold"))

#--------Temporal: Plotting Top 10 Most Abundant by Day of Life---------------------------
top10_each %>%
  ggplot(aes(fill = Genus,
             y = Genus_Count, 
             x = Day_of_Life)) +                          #This shows a temporal view seperated by time
  geom_bar(position = "fill", 
           stat = "identity") +
  labs(title = "Total Genus Abundance per Sample",
       x = "Date",
       y = "Frequency of Abundance",
       fill = "Genus") +
  scale_y_continuous(labels = percent) +
  annotate("text", x = 23, y = -0.038, label= "9/11/2019", size = 4, angle = 24, hjust = .6) + 
  annotate("text", x = 35, y = -0.038, label= "9/23/2019", size = 4, angle = 24, hjust = .6) +
  annotate("text", x = 38, y = -0.038, label= "9/26/2019", size = 4, angle = 24, hjust = .6) +
  annotate("text", x = 42, y = -0.038, label= "9/30/2019", size = 4, angle = 24, hjust = .6) +
  annotate("text", x = 48, y = -0.038, label= "10/06/2019", size = 4, angle = 24, hjust = .6) +
  theme(plot.title = element_text(size = 18, face = "bold.italic"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
