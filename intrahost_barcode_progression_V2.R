library(tidyverse)

#set working directory and read in files:
#change path to location of files.
#sample files
getwd()
setwd("/Users/emilyfitzmeyer/Desktop/BC_WNV/ya_bcWNV/all_counts/")
column_classes <- c("c", "n")

#make sure you are reading in samples from the same dpi and mosquito!
ser <- read_delim("Ifnar_Ser_F80-C7_true_barcodes.txt", 
                 col_types = column_classes, col_names = c("barcode", "count"))
spl <- read_delim("Ifnar_Spl_F80-F7_true_barcodes.txt", 
                 col_types = column_classes, col_names = c("barcode", "count"))
ute <- read_delim("IFNAR_Ute_80_true_barcodes.txt", 
                 col_types = column_classes, col_names = c("barcode", "count"))
dec <- read_delim("Ifnar_Dec_F80-3-F8_true_barcodes.txt", 
                 col_types = column_classes, col_names = c("barcode", "count"))
pla <- read_delim("Ifnar_Pla_F80-3-G9_true_barcodes.txt", 
                 col_types = column_classes, col_names = c("barcode", "count"))

names <- c("ser","spl", "ute", "dec", "pla")

list1 <- list(ser, spl, ute, dec, pla)

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
  select(1, 3, 5, 7, 9, 11)

#arrange by MG or input 
graph_data <- graph_data %>%
  arrange(desc(freq_spl))

#generate separate dataframes so each element of the X axis gets its own rank

sepSer <- graph_data[c("freq_ser")]
sepSer <- sepSer %>%
  mutate(sample = 0) %>%
  mutate(barcode_id = 1:nrow(sepSer)) %>%
  dplyr::rename("freq" = "freq_ser")

sepSpl <- graph_data[c("freq_spl")]
sepSpl <- sepSpl %>%
  mutate(sample = 1) %>%
  mutate(barcode_id = 1:nrow(sepSpl)) %>%
  dplyr::rename("freq" = "freq_spl")

sepUte <- graph_data[c("freq_ute")]
sepUte <- sepUte %>%
  mutate(sample = 2) %>%
  mutate(barcode_id = 1:nrow(sepUte)) %>%
  dplyr::rename("freq" = "freq_ute")

sepDec <- graph_data[c("freq_dec")]
sepDec <- sepDec %>%
  mutate(sample = 3) %>%
  mutate(barcode_id = 1:nrow(sepDec)) %>%
  dplyr::rename("freq" = "freq_dec")

sepPla <- graph_data[c("freq_pla")]
sepPla <- sepPla %>%
  mutate(sample = 4) %>%
  mutate(barcode_id = 1:nrow(sepPla)) %>%
  dplyr::rename("freq" = "freq_pla")

#bind these dataframes and arrange desc by barcode_id
frequencies <- rbind(sepSer,sepSpl, sepUte, sepDec, sepPla) %>%
  arrange(desc(barcode_id)) %>%
  mutate(barcode_id = as.factor(barcode_id)) #%>%
  #mutate(volume = as.factor(sample))

#graph
x<-(rep(c("red","green","blue","yellow","purple","grey","orange","brown","turquoise","pink","cyan","darkgrey"),1500))

breaks <- c("Ser", "Spl", "Ute", "Dec", "Pla")

ggplot(frequencies, aes(x = sample, y = freq, fill = barcode_id)) + 
  geom_area(alpha = 0.6, show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(panel.grid.major.y = element_line(colour = "light grey"), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(axis.title.x = element_text(size = 1), axis.title.y = element_text(size = 1)) +
  theme(plot.margin = unit(c(0.3,1,0.3,0.3), "cm")) +
  ggtitle("IFNAR_80.3") +
  theme(plot.title = element_text(vjust = 4)) +
  labs(x = "", y = "") +
  scale_x_continuous(breaks = c(0,1,2,3,4), labels = stringr::str_wrap(breaks, width = 8), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1, by=0.2), expand = c(0,0)) +
  scale_fill_manual(values = x)

#setwd("/Users/emilyfitzmeyer/Desktop/prog_plots/Supplemental/")
#ggsave("cxt_9_12dpi.png", plot = last_plot(), device = png(), scale = 1, width = 6.5, height = 5, dpi = 300)

