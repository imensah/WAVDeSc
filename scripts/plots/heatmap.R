#libraries
library(ggplot2)
library(reshape2)
library(viridis)

# Load NONUMI_12 data
data <- read.csv("/home/irenekay/Documents/thesis_scrna/nonumi_12.csv")

# Melt the data for use with ggplot
melted_data <- melt(data, id.vars = 'MODEL')

# Create the heatmap with a gradient from cyan to magenta
ggplot(melted_data, aes(x = MODEL, y = variable, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = value), color = "white", size = 9, vjust = 1, fontface = "bold") +
  scale_fill_gradient(low = "lightblue", high = "maroon", name = "Value") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 20),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold")) +
  labs(title = "")


# LOAD DATA: UMI_12
data <- read.csv("/home/irenekay/Documents/thesis_scrna/umi_12.csv")

# Melt the data for use with ggplot
melted_data <- melt(data, id.vars = 'MODEL')

# Create the heatmap with a gradient from cyan to magenta
ggplot(melted_data, aes(x = MODEL, y = variable, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = value), color = "white", size = 9, vjust = 1, fontface = "bold") +
  scale_fill_gradient(low = "lightblue", high = "maroon", name = "Value") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 20),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold")) +
  labs(title = "")

# LOAD DATA: NON-UMI_34
data <- read.csv("/home/irenekay/Documents/thesis_scrna/nonumi_34.csv")

# Melt the data for use with ggplot
melted_data <- melt(data, id.vars = 'MODEL')

# Create the heatmap with a gradient from cyan to magenta
ggplot(melted_data, aes(x = MODEL, y = variable, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = value), color = "white", size = 9, vjust = 1, fontface = "bold") +
  scale_fill_gradient(low = "lightblue", high = "maroon", name = "Value") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 20),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold")) +
  labs(title = "")

# LOAD DATA: UMI_34
data <- read.csv("/home/irenekay/Documents/thesis_scrna/umi_34.csv")

# Melt the data for use with ggplot
melted_data <- melt(data, id.vars = 'MODEL')

# Create the heatmap with a gradient from cyan to magenta
ggplot(melted_data, aes(x = MODEL, y = variable, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = value), color = "white", size = 9, vjust = 1, fontface = "bold") +
  scale_fill_gradient(low = "lightblue", high = "maroon", name = "Value") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 20),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold")) +
  labs(title = "")




