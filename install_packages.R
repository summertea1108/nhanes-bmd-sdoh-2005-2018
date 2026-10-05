# Run ONCE on an OSC login node (compute nodes may not have internet).
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
pkgs <- c("haven", "dplyr", "tidyr", "purrr", "readr", "survey", "mitools",
          "mice", "ranger", "iml", "ggplot2")
install.packages(setdiff(pkgs, rownames(installed.packages())),
                 lib = lib, repos = "https://cloud.r-project.org")
cat("Missing after install:", setdiff(pkgs, rownames(installed.packages())), "\n")
