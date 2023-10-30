library(tidyverse)

#set working directory and read in files:
#change path to location of files.
#sample files
getwd()
setwd("/Users/emilyfitzmeyer/Desktop/BC_WNV/ya_bcWNV/all_counts/")
column_classes <- c("c", "n")

#read in all Ute, or Dec, or Pla samples
spl1 <- read_delim("Mavs_Spl_F71-E7_true_barcodes.txt", 
                   col_types = column_classes, col_names = c("barcode", "count"))
ser1 <- read_delim("Mavs_Ser_F71-B7_true_barcodes.txt", 
                   col_types = column_classes, col_names = c("barcode", "count"))
ute1 <- read_delim("Mavs_Ute_F71-G7_true_barcodes.txt", 
                 col_types = column_classes, col_names = c("barcode", "count"))
#dec1 <- read_delim("WT_Dec_F260-1-G2_true_barcodes.txt", 
#                 col_types = column_classes, col_names = c("barcode", "count"))
#dec2 <- read_delim("Mavs_Dec_F68-2-H7_true_barcodes.txt", 
#                 col_types = column_classes, col_names = c("barcode", "count"))
#dec3 <- read_delim("Mavs_Dec_F68-3-A8_true_barcodes.txt", 
#                  col_types = column_classes, col_names = c("barcode", "count"))
dec4 <- read_delim("Mavs_Dec_F71-4-C8_true_barcodes.txt", 
                  col_types = column_classes, col_names = c("barcode", "count"))
dec5 <- read_delim("MAVS_Dec2_71.5_true_barcodes.txt", 
                  col_types = column_classes, col_names = c("barcode", "count"))
dec6 <- read_delim("Mavs_Dec_F71-6-D8_true_barcodes.txt", 
                  col_types = column_classes, col_names = c("barcode", "count"))
#dec7 <- read_delim("WT_Dec_F262-7-C4_true_barcodes.txt", 
#                  col_types = column_classes, col_names = c("barcode", "count"))
#dec8 <- read_delim("WT_Dec_F262-8-D4_true_barcodes.txt", 
#                 col_types = column_classes, col_names = c("barcode", "count"))
#pla1 <- read_delim("WT_Dec_F260-1-G2_true_barcodes.txt", 
#                   col_types = column_classes, col_names = c("barcode", "count"))
#pla2 <- read_delim("Mavs_Pla_F68-2-A9_true_barcodes.txt", 
#                   col_types = column_classes, col_names = c("barcode", "count"))
#pla3 <- read_delim("Mavs_Pla_F68-3-B9_true_barcodes.txt", 
#                   col_types = column_classes, col_names = c("barcode", "count"))
pla4 <- read_delim("Mavs_Pla_F71-4-D9_true_barcodes.txt", 
                   col_types = column_classes, col_names = c("barcode", "count"))
pla5 <- read_delim("MAVS_Pla2_71.5_true_barcodes.txt", 
                   col_types = column_classes, col_names = c("barcode", "count"))
pla6 <- read_delim("Mavs_Pla_F71-6-E9_true_barcodes.txt", 
                   col_types = column_classes, col_names = c("barcode", "count"))
#pla7 <- read_delim("WT_Pla_F262-7-G6_true_barcodes.txt", 
#                   col_types = column_classes, col_names = c("barcode", "count"))
#pla8 <- read_delim("WT_Dec_F262-8-D4_true_barcodes.txt", 
#                   col_types = column_classes, col_names = c("barcode", "count"))

names <- c("spl1", "ser1", "ute1", "dec4", "dec5", "dec6",
           "pla4", "pla5", "pla6")

#, "dec1", "dec2", "dec3", "dec4", "dec5", "dec6", "dec7", "dec8"
#, "pla1", "pla2", "pla3", "pla4", "pla5", "pla6", "pla7", "pla8"

list1 <- list(spl1, ser1, ute1, dec4, dec5, dec6,
              pla4, pla5, pla6)

#, dec1, dec2, dec3, dec5, dec6, dec7, dec8
#, pla1, pla2, pla3, pla5, pla6, pla7, pla8

names(list1) <- names

#create function to generate frequency columns for each df and preserve + append df name to numeric columns
barcode_frequency <- function(x, name) {
  totalumis <- sum(x[["count"]])
  x[["freq"]] <- ((x[["count"]])/totalumis)
  names(x) <- paste(names(x), name, sep = "_")
  colnames(x)[1] = "barcode"
  print(x)
}

#map function over list
list2 <- map2(list1, names(list1), barcode_frequency)

#merge dfs by barcode
all_samples <- list2 %>%
  reduce(full_join, by = "barcode")

#make NA entries '0'
all_samples[is.na(all_samples)]<-0

#select frequency columns for graph data
graph_data <- all_samples %>%
  select(1, 3, 5, 7, 9, 11, 13, 15, 17, 19) 
#1:39 = 18 total samples

#arrange desc by input
graph_data <- graph_data %>%
  arrange(desc(freq_spl1))

#generate separate dataframes so each element of the X axis gets its own rank
#remember to update sample = # if commenting out samples!!!!!!!
sepSpl1 <- graph_data[c("freq_spl1")]
sepSpl1 <- sepSpl1 %>%
  mutate(sample = 0) %>%
  mutate(barcode_id = 1:nrow(sepSpl1)) %>%
  dplyr::rename("freq" = "freq_spl1")

sepSer1 <- graph_data[c("freq_ser1")]
sepSer1 <- sepSer1 %>%
  mutate(sample = 1) %>%
  mutate(barcode_id = 1:nrow(sepSer1)) %>%
  dplyr::rename("freq" = "freq_ser1")

sepUte1 <- graph_data[c("freq_ute1")]
sepUte1 <- sepUte1 %>%
  mutate(sample = 2) %>%
  mutate(barcode_id = 1:nrow(sepUte1)) %>%
  dplyr::rename("freq" = "freq_ute1")

#sepDec1 <- graph_data[c("freq_dec1")]
#sepDec1 <- sepDec1 %>%
#  mutate(sample = 3) %>%
#  mutate(barcode_id = 1:nrow(sepDec1)) %>%
#  dplyr::rename("freq" = "freq_dec1")

#sepDec2 <- graph_data[c("freq_dec2")]
#sepDec2 <- sepDec2 %>%
#  mutate(sample = 3) %>%
#  mutate(barcode_id = 1:nrow(sepDec2)) %>%
#  dplyr::rename("freq" = "freq_dec2")

#sepDec3 <- graph_data[c("freq_dec3")]
#sepDec3 <- sepDec3 %>%
#  mutate(sample = 4) %>%
#  mutate(barcode_id = 1:nrow(sepDec3)) %>%
#  dplyr::rename("freq" = "freq_dec3")

sepDec4 <- graph_data[c("freq_dec4")]
sepDec4 <- sepDec4 %>%
  mutate(sample = 3) %>%
  mutate(barcode_id = 1:nrow(sepDec4)) %>%
  dplyr::rename("freq" = "freq_dec4")

sepDec5 <- graph_data[c("freq_dec5")]
sepDec5 <- sepDec5 %>%
  mutate(sample = 4) %>%
  mutate(barcode_id = 1:nrow(sepDec5)) %>%
  dplyr::rename("freq" = "freq_dec5")

sepDec6 <- graph_data[c("freq_dec6")]
sepDec6 <- sepDec6 %>%
  mutate(sample = 5) %>%
  mutate(barcode_id = 1:nrow(sepDec6)) %>%
  dplyr::rename("freq" = "freq_dec6")

#sepDec7 <- graph_data[c("freq_dec7")]
#sepDec7 <- sepDec7 %>%
#  mutate(sample = 8) %>%
#  mutate(barcode_id = 1:nrow(sepDec7)) %>%
#  dplyr::rename("freq" = "freq_dec7")

#sepDec8 <- graph_data[c("freq_dec8")]
#sepDec8 <- sepDec8 %>%
#  mutate(sample = 9) %>%
#  mutate(barcode_id = 1:nrow(sepDec8)) %>%
#  dplyr::rename("freq" = "freq_dec8")

#sepPla1 <- graph_data[c("freq_pla1")]
#sepPla1 <- sepPla1 %>%
#  mutate(sample = 9) %>%
#  mutate(barcode_id = 1:nrow(sepPla1)) %>%
#  dplyr::rename("freq" = "freq_pla1")

#sepPla2 <- graph_data[c("freq_pla2")]
#sepPla2 <- sepPla2 %>%
#  mutate(sample = 7) %>%
#  mutate(barcode_id = 1:nrow(sepPla2)) %>%
#  dplyr::rename("freq" = "freq_pla2")

#sepPla3 <- graph_data[c("freq_pla3")]
#sepPla3 <- sepPla3 %>%
#  mutate(sample = 8) %>%
#  mutate(barcode_id = 1:nrow(sepPla3)) %>%
#  dplyr::rename("freq" = "freq_pla3")

sepPla4 <- graph_data[c("freq_pla4")]
sepPla4 <- sepPla4 %>%
  mutate(sample = 6) %>%
  mutate(barcode_id = 1:nrow(sepPla4)) %>%
  dplyr::rename("freq" = "freq_pla4")

sepPla5 <- graph_data[c("freq_pla5")]
sepPla5 <- sepPla5 %>%
  mutate(sample = 7) %>%
  mutate(barcode_id = 1:nrow(sepPla5)) %>%
  dplyr::rename("freq" = "freq_pla5")

sepPla6 <- graph_data[c("freq_pla6")]
sepPla6 <- sepPla6 %>%
  mutate(sample = 8) %>%
  mutate(barcode_id = 1:nrow(sepPla6)) %>%
  dplyr::rename("freq" = "freq_pla6")

#sepPla7 <- graph_data[c("freq_pla7")]
#sepPla7 <- sepPla7 %>%
#  mutate(sample = 16) %>%
#  mutate(barcode_id = 1:nrow(sepPla7)) %>%
#  dplyr::rename("freq" = "freq_pla7")

#sepPla8 <- graph_data[c("freq_pla8")]
#sepPla8 <- sepPla8 %>%
#  mutate(sample = 17) %>%
#  mutate(barcode_id = 1:nrow(sepPla8)) %>%
#  dplyr::rename("freq" = "freq_pla8")

#bind these dataframes and arrange desc by barcode_id
frequencies <- rbind(sepSpl1, sepSer1, sepUte1, sepDec4, sepDec5, sepDec6,
                     sepPla4, sepPla5, sepPla6) %>%
  arrange(desc(barcode_id)) %>%
  mutate(barcode_id = as.factor(barcode_id))
#, sepDec1, sepDec2, sepDec3, sepDec4, sepDec5, sepDec6, sepDec7, sepDec8
#, sepPla1, sepPla2, sepPla3, sepPla4, sepPla5, sepPla6, sepPla7, sepPla8

#graph - note if getting error "insufficient values in manual scale" increase the value of the
#number after the color vector. 
x<-(rep(c("red","green","blue","yellow","purple","grey","orange","brown","turquoise","pink",
          "cyan","darkgrey"),10000))

breaks <- c("Spl", "Ser", "Ute", "Dec4", "Dec5", "Dec6",
            "Pla4", "Pla5", "Pla6")
#"Dec1", "Dec2", "Dec3", "Dec7", "Dec8"
#"Pla1", "Pla2", "Pla3", "Pla7", "Pla8"

#have to include 'stat="identity"' when plotting this data with geom_bar
ggplot(frequencies, aes(x = sample, y = freq, fill = barcode_id)) + 
  geom_bar(stat = "identity", alpha = 0.6, show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        panel.background = element_blank()) +
  theme(panel.grid.major.y = element_line(colour = "light grey"), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(axis.title.x = element_text(size = 1), axis.title.y = element_text(size = 1)) +
  theme(plot.margin = unit(c(0.5,0.1,0.1,0.1), "cm")) +
  ggtitle("MAVS_71") +
  theme(plot.title = element_text(vjust = 4)) +
  labs(x = "", y = "") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8), labels = breaks, expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1, by=0.2), expand = c(0,0)) +
  scale_fill_manual(values = x)

setwd("/Users/emilyfitzmeyer/Desktop/")
ggsave("MAVS_71.png", plot = last_plot(), device = png(), scale = 1, width = 12, height = 5, dpi = 300)

