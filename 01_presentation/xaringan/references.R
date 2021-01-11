if (!is.element("RefManageR", installed.packages()[,1])) install.packages("RefManageR")
if (!is.element("tidyverse", installed.packages()[,1])) install.packages("tidyverse")

# References function --------------------------------------------------
library(RefManageR)
library(tidyverse)

references <- function(bib_file, select = FALSE, keys = NULL){
  bib <- ReadBib(bib_file, check = FALSE)
  
  if(select){
    print(bib[key = keys], 
          .opts = list(check.entries = FALSE, 
                       style = "html", 
                       bib.style = "authoryear"))
  } else {
    print(bib, .opts = list(check.entries = FALSE, 
                            style = "html", 
                            bib.style = "authoryear"))
    
  }
  
}