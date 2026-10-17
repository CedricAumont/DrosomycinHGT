# CMH test within order

data_array1 <- array(c(
  72, 0, 32, 12,
  71, 1, 35, 9,
  65, 7, 33, 11,
  71, 6, 38, 6,
  49, 23, 30, 15,
  66, 7, 40, 4,
  61, 12, 40, 6,
  71, 1, 36, 8
),
dim = c(2,2,8),
dimnames = list(
  Response = c("id","non-id"),
  Condition = c("IPT","Drs"),
  Comparisons = c("K","M","B","Br","P","E","C","D")
)
)
cmh_test1 = mantelhaen.test(data_array1, alternative =  "less")
cmh_test1

# 	Mantel-Haenszel chi-squared test with continuity correction

# data:  data_array1
# Mantel-Haenszel X-squared = 19.104, df = 1, p-value = 1
# alternative hypothesis: true common odds ratio is less than 1
# 95 percent confidence interval:
#   0.000000 3.216705
# sample estimates:
#   common odds ratio 
# 2.328818 

# CMH test for order consensus pairwise comparison

data_array <- array(c(
  58,14,36,8,
  57,15,36,8,
  49,23,38,6,
  49,23,32,12,
  53,19,38,6,
  50,22,38,6,
  46,26,38,6,
  46,26,30,14,
  62,10,37,7,
  52,20,41,3,
  47,25,34,10,
  53,19,41,3,
  51,21,40,4,
  50,22,41,3,
  51,21,29,15,
  51,21,40,4,
  47,25,36,8,
  51,21,37,7,
  47,25,39,5,
  49,23,40,4,
  49,23,32,12,
  48,24,35,9,
  46,26,39,5,
  46,26,41,3,
  45,27,44,0,
  46,26,32,12,
  40,32,34,10,
  44,28,34,10,
  43,29,35,9,
  43,29,29,15,
  55,17,42,2,
  47,25,39,5,
  46,26,29,15,
  45,27,41,3,
  45,27,31,13,
  69,3,32,12
),
dim = c(2,2,36),
dimnames = list(
  Response = c("id","non-id"),
  Condition = c("RDH","Drs"),
  Comparisons = c("KM", "KR", "KB", "KP", "KE", "KC",
                  "KD", "Kbu", "MR", "MB", "MP", "ME", "MC", 
                  "MD", "Mbu", "RB", "RP", "RE", "RC", "RD", "Rbu",
                  "BP", "BE", "BC", "BD", "Bbu", "PE", "PC", "PD",
                  "Pbu", "EC", "ED", "Ebu", "CD", "Cbu", "Dbu"
)
)
)
cmh_test = mantelhaen.test(data_array, alternative = "less")
cmh_test
# 	Mantel-Haenszel chi-squared test with continuity correction
# 
# data:  data_array
# Mantel-Haenszel X-squared = 102.97, df = 1, p-value < 2.2e-16
# alternative hypothesis: true common odds ratio is less than 1
# 95 percent confidence interval:
#   0.0000000 0.5168535
# sample estimates:
#   common odds ratio 
# 0.4538272 
