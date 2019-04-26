
install.packages('fmsb')
library(fmsb)

M <- matrix(c(759, 360, 518, 363), ncol = 2)
oddsratio(M)

# The number of DNA regions in each category = 6798

#                                    | The number of triplexes predicted |
#                                    | by Triplexator is statistically   |
#                                    | significant                       |
#                                    |-----------------------------------|
#                                    |     YES        |       NO         |
#                                    |----------------|------------------|
# The MEG3-DNA interaction was | YES |     3825       |       2973       |
# identified by ChOP-seq       | NO  |     617        |       6181       |
# 
M <- matrix(
  c(3825, 6798-3825,
    617,  6798-617),
  ncol = 2, byrow = TRUE)
oddsratio(M)
