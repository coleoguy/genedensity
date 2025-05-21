# data
# indices are divergence scores from 0 to 16
# values are # of nucleotides in each bin
vec <- c(3, 10, 140, 171, 231, 254, 213, 130, 64, 22, 6, 1, 0, 0, 0, 0, 0)
plot(vec) # plot for reference

# current method
which(cumsum(vec) > sum(vec)/2)[1]

# using the median's divergence as repeat age
which(vec == median(vec))

# using the divergence of the number closest to the mean as repeat age
res <- abs(mean(vec) - vec)
which(min(res) == res)

# expand histogram data into individual observations
idx <- 1:length(vec)-1
values <- rep(idx, vec)

mean(values) # find mean
median(values) # find median


