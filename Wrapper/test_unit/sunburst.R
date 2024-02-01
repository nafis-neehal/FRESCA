# Install and load required packages
library(plotly)

# Create a data frame with hierarchical data
# id will be replaced by labels ****
# parent should be the values in id ****


data <- data.frame(
  id =    c( "Male", "Female", "Male/40-59", "Male/59+",  "Female/40-59",  "Female/59+",  "Male/40-59/Hispanic",  "Female/59+/Others"),         
  label = c( "Male", "Female", "40-59",      "59+",       "40-59",         "59+",         "Hispanic",             "Others"),
  parent = c("",     "",       "Male",       "Male",      "Female",        "Female",      "Male/40-59",           "Female/59+"),                
  value = c(50, 50, 25, 25, 25, 25, 5, 5)
)

# Create a sunburst plot
plot <- plot_ly(
  data,
  ids = ~id,
  labels = ~label,
  parents = ~parent,
  values = ~value,
  type = "sunburst",
  branchvalues="total"
)

# Customize the sunburst plot
plot <- plot %>% layout(
  sunburstcolorway = c("#636efa", "#EF553B", "#00cc96", "#ab63fa", "#19d3f3", "#FFA15A", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"),
  margin = list(t = 0, l = 0, r = 0, b = 0),
  sunburst = list(
    branchvalues = "total",
    labels = list(
      rotate = 0,
      font = list(size = 14)
    )
  )
)

# Display the sunburst plot
plot





#,
#"Female/40-59/NH Asian", "Female/40-59/NH Black", "Female/40-59/NH White", "Female/40-59/Hispanic", "Female/40-59/Others",
#"Female/59+/NH Asian", "Female/59+/NH Black", "Female/59+/NH White", "Female/59+/Hispanic", "Female/59+/Others"

#,
#"40-59","40-59","40-59","40-59","40-59","59+","59+","59+","59+","59+"
#,
#"NH Asian", "NH Black", "NH White", "Hispanic", "Others","NH Asian", "NH Black", "NH White", "Hispanic", "Others"
#, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10

#"Male/40-59/NH Asian", "Male/40-59/NH Black", "Male/40-59/NH White", "Male/40-59/Hispanic", "Male/40-59/Others",
#"Male/59+/NH Asian", "Male/59+/NH Black", "Male/59+/NH White", "Male/59+/Hispanic", "Male/59+/Others"
#, "40-59",    "40-59",    "40-59",    "40-59",    "40-59",  "59+",      "59+",      "59+",      "59+",      "59+
#"NH Asian", "NH Black", "NH White", "Hispanic", "Others", "NH Asian", "NH Black", "NH White", "Hispanic", "Others"


# # Create a data frame with hierarchical data
# data <- data.frame(
#   #ids are for each labels so that labels are not repeated
#   id = c("Male", "Female", "Male/40-59", "Male/59+", "Female/40-59", "Female/59+", "Male/40-59/H"), 
#   parent = c("",     "",       "Male",  "Male", "Female", "Female", "40-59"),
#   label =  c("Male", "Female", "40-59", "59+",  "40-59",  "59+",    "Hispanic"),
#   value = c(50, 50, 25, 25, 25, 25, 5)
# )


# # Create a sunburst plot
# plot <- plot_ly(
#   data,
#   ids = ~id,
#   labels = ~label,
#   parents = ~parent,
#   values = ~value,
#   type = "sunburst",
#   branchvalues = 'total'
# )
# 
# # Display the sunburst plot
# plot
