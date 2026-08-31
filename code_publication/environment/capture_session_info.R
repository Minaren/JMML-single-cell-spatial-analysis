# Run from the code_publication root after installing the required packages.
dir.create("output", showWarnings = FALSE, recursive = TRUE)
sink(file.path("output", "sessionInfo.txt"))
print(sessionInfo())
sink()
source(file.path("environment", "packages.R"))
