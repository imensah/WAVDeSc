

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggrepel)

# Function to plot and customize a line graph

plot_line_graph <- function(data) {
  plot <- ggplot(data, aes(x = variable, y = value, group = Model)) +
    geom_line(aes(color = Model), size = 1.5) + # Increase line width
    scale_y_continuous(limits = c(0.5, 1)) +   # Set y-axis limits
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
          panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
          panel.grid.major = element_line(size = 0.1, color = "gray"),  # Remove major gridlines
          panel.grid.minor = element_line(size = 0.1, color = "gray"),  # Remove minor gridlines
          axis.text.x.bottom = element_text(size = 28, face = "bold"),
          axis.text.y = element_text(size = 25, face = "bold"))
  
  data_ends <- data %>% filter(variable == "F1.Score")
  
  plot + geom_text_repel(
    aes(label = Model), data = data_ends,
    color = "#c00000", size = 5, direction = "y", fontface = "bold", xlim = c(4, NA),
    hjust = 0,
    segment.size = .7,
    segment.alpha = -0.1,
    segment.curvature = 0.5,
    segment.ncp = 3
  )
}

# Load and process the first dataset umi12
Data1 <- read.csv('/home/irenekay/Downloads/IRENEKAY/Isabel/line_graph/umi/pop12.csv')
#Data1[4, 1] <- 'WAVDeSc'
#remove specificity
dd = Data1[,-4]

#remove FDR
dr = dd[,-6]

#remove db6
#dc = dr[-5,]

data1 <- melt(dr)


# Plot the first graph
plot_line_graph(data1)




# Load and process the second dataset umi34
Data2 <- read.csv('/home/irenekay/Downloads/IRENEKAY/Isabel/line_graph/umi/pop34.csv')
#Data2[4, 1] <- 'WAVDeSc'
#remove specificity
dd = Data2[,-4]
dd
#remove FDR
dr = dd[,-6]
dr
#remove db6
#dc = dr[-5,]
data2 <- melt(dr)

# Plot the second graph
plot_line_graph(data2)

# Load and process the third dataset nonumi12
Data3 = read.csv('/home/irenekay/Downloads/IRENEKAY/Isabel/line_graph/nonumi/npop12.csv')
#Data3[4,1] = 'WAVDeSc'
#remove specificity
dd = Data3[,-4]
dd
#remove FDR
dr = dd[,-6]
dr
#remove db6
#dc = dr[-5,]
data3 = melt(dr)

# Plot the third graph
plot_line_graph(data3)

# Load and process the fourth dataset nonumi34
Data4 = read.csv('/home/irenekay/Downloads/IRENEKAY/Isabel/line_graph/nonumi/npop34.csv')
#Data4[4,1] = 'WAVDeSc'
#remove specificity
dd = Data4[,-4]
dd
#remove FDR
dr = dd[,-6]
dr
#remove db6
#dc = dr[-5,]
data4 = melt(dr)

# Plot the third graph
plot_line_graph(data4)
