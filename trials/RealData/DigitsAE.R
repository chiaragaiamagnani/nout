# DIGITS PRE-PROCESSING

library(reshape2)
library(ggplot2)
library(FactoMineR)

pen.tra = read.table("~/nout/trials/RealData/Datasets/Dataset digits/pendigits.tra", sep = ",")
pen.tes = read.table("~/nout/trials/RealData/Datasets/Dataset digits/pendigits.tes", sep = ",")

pen = rbind(pen.tra, pen.tes)

names(pen) = c(paste0(c("X", "Y"), rep(1:8, each = 2)), "digit")
class(pen$digit)

pen$digit = factor(pen$digit)
class(pen$digit)

table(pen$digit)
ggplot(pen, aes(digit, fill = digit)) +
  geom_bar() +
  guides(fill = FALSE)

tab = apply(pen[-17], 2, function(x) {
  return (c(Mean = mean(x),
            Median = median(x),
            StdDev = sd(x)))
})
round(t(tab), 2)

# Boxplot of each coordinate
for (j in 1:16)
  boxplot(pen[,j])
pen.X = pen[seq(1, 15, 2)]
names(pen.X)
pen.X.melt = melt(pen.X)
knitr::kable(head(pen.X.melt))
ggplot(pen.X.melt, aes(variable, value)) +
  geom_boxplot()

pen.Y = pen[seq(2, 16, 2)]
names(pen.Y)
pen.Y.melt = melt(pen.Y)
knitr::kable(head(pen.Y.melt))
ggplot(pen.Y.melt, aes(variable, value)) + geom_boxplot()


# Now, we have the basis of a function that can reproduc each digit. We want the following constraints for this one :

# Parameters:
#   - v: a vector with the 16 values of each point (if it is a data.frame, transform it as simple vector)
#   - n: a numeric value (with NULL value by default), denoting the number to drawn (possibly not pass)
#   - point: a boolean vector (with FALSE value by default), indicating if the point number has to be added to the graph
#   - add: a boolean vector (with FALSE value by default), indicating if the plot is to added to the previous one
#   - color: a specific color (with "black" value by default)
# Output:
#   - line representing the drawn
#   - number drawn as title (if known)
#   - point number (if wanted)

drawn <- function(v, n = NULL, point = FALSE, add = FALSE, color = "black") {
  # Transformation of the data.frame if needed
  if (is.data.frame(v))
    v = unlist(v)
  # extract x and y coordinates
  x = v[seq(1, 15, by = 2)]
  y = v[seq(2, 16, by = 2)]
  # optimize space into graphics in reducing margin (sse ?par for more information)
  opar = par(mar = c(0, 0, 2, 0) + .1)
  if (!add) { # Create a graphic
    plot(x, y,
         # Specify limits is a way to have always the same frame for plotting lines
         xlim = c(-5, 105), ylim = c(-5, 105),
         # Do not show axes
         axes = FALSE,
         # If point is TRUE, we add a space (with pch = " ") at each point
         # If not, draw a line
         type = ifelse(point, "b", "l"),
         pch = " ",
         # Specify color (black by default)
         col = color,
         # Add a title (NULL by default)
         main = n)
    if (point) text(x, y, 1:8)
  } else { # Add line to the plot
    # lines() add lines to an existing plot (the last produce)
    lines(x, y,
          # same comment than before
          xlim = c(-5, 105), ylim = c(-5, 105),
          type = ifelse(point, "b", "l"),
          pch = " ",
          col = color)
  }
  par(opar)
}



first.0 = subset(pen, digit == 0)[1,]
drawn(first.0)
drawn(first.0, n=0)
drawn(first.0, point=TRUE)

first.1 = subset(pen, digit == 1)[1,]
drawn(first.1, add = TRUE, col = "red")

par(mfrow = c(2, 5))
for (i in 0:9) {
  s = subset(pen, digit == i)[1,]
  drawn(s[-16], n = i)
}


# We now want to represent the average digit,
# i.e. a line produce by the average value for each coordinate.
# To see if there are some area that are not use for a particular digit,
# we represent also all the drawn in a very light gray.

# First, we present the mean value of each coordinate for each digit.
# We can do this step with apply() (and tapply()) or with aggregate().

apply(pen[-17], 2, tapply, pen$digit, mean)
agg = aggregate(. ~ digit, pen, mean)
knitr::kable(agg, digits = 2)

# We now represent the average digit.
# To see if there are some area that are not use for a particular digit,
# we represent also all the drawn in a very light gray.

par(mfrow = c(2, 5))
for (i in 0:9) {
  s = subset(pen, digit == i)
  # Computing the average values
  mean = apply(s[-17], 2, mean)
  # Drawn the first drawn
  drawn(s[1,1:16], col = "gray90", n = i)
  # Add all other drawn
  for (j in 2:nrow(s))
    drawn(s[j,1:16], col = "gray90", add = TRUE)
  # Last, add the average line
  drawn(mean[-17], add = TRUE)
}


# PCA
res = PCA(pen, quali.sup = 17, graph = FALSE)
eig = data.frame(comp = 1:16,
                 setNames(res$eig, c("eigenvalue", "percent", "cum.percent")))
ggplot(eig) +
  geom_bar(aes(comp, percent), stat = "identity") +
  stat_summary(aes(comp, cum.percent, group = 1),
               fun.y = sum, geom = "line")
res2 = data.frame(res$ind$coord, digit = pen$digit)
ggplot(res2, aes(Dim.1, Dim.2, color = digit)) +
  geom_point()


# GGLPOT
# We can also represent digit drawn with the ggplot() function.
# We have to use the melt() function, and to be careful on how we do this.

# For the first 0

f0 = data.frame(x = unlist(first.0[seq(1, 15, 2)]),
                y = unlist(first.0[seq(2, 16, 2)]))
ggplot(f0, aes(x, y)) +
  geom_path() +
  theme_void()


# For the average digit, and all drawn.

# Construct the data for all drawn
pen$id = 1:nrow(pen)
pen.melt = melt(pen, id.vars = c("id", "digit"))
pen.melt$coord = substr(pen.melt$variable, 1, 1)
pen.melt$point = substr(pen.melt$variable, 2, 2)
pen.melt.X = subset(pen.melt, coord == "X", -c(variable, coord))
names(pen.melt.X) = c("id", "digit", "x", "point")
pen.melt.Y = subset(pen.melt, coord == "Y", -c(variable, coord))
names(pen.melt.Y) = c("id", "digit", "y", "point")
pen.merge = merge(pen.melt.X, pen.melt.Y)

# Construct the data for average

agg.melt = melt(agg, id.vars = "digit")
agg.melt$coord = substr(agg.melt$variable, 1, 1)
agg.melt$point = substr(agg.melt$variable, 2, 2)
agg.melt.X = subset(agg.melt, coord == "X", -c(variable, coord))
names(agg.melt.X) = c("digit", "x", "point")
agg.melt.Y = subset(agg.melt, coord == "Y", -c(variable, coord))
names(agg.melt.Y) = c("digit", "y", "point")
agg.merge = merge(agg.melt.X, agg.melt.Y)

# Plot it

ggplot(pen.merge, aes(x, y)) +
  geom_path(aes(group = id), color = "gray90") +
  facet_wrap(~ digit) +
  geom_path(data = agg.merge) +
  theme_void()




data = readMat("~/nout/trials/RealData/Datasets/Dataset digits/pendigits.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

length(out_ind)


