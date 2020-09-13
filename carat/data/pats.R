set.seed(201911)
pats <- data.frame("gender" = 
                     sample(c("Female", "Male"), 1000, TRUE), 
                   "employment status" = 
                     sample(c("unemp", "part.", "full."), 1000, TRUE), 
                   "income" = 
                     sample(c(">= 1w", "<= 0.5w", "0.5~1w"), 1000, TRUE), 
                   "marriage status" = 
                     sample(c("unmarried", "married", "divorced"), 1000, TRUE), 
                   stringsAsFactors = TRUE)
