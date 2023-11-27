setwd("~/repos/official_repos/mappoly2/")
# Get a list of all files
all_files <- list.files(path = ".", recursive = TRUE)

# Filter for R and C++ source files
source_files <- all_files[grep("\\.R$|\\.cpp$", all_files)]

# Initialize empty list to store function names
function_names <- list()

# Define regex patterns for R and C++ function declarations
r_function_pattern <- "^.*<-\\s*function"
cpp_function_pattern <- "^[a-zA-Z_]+\\s+[a-zA-Z_]+\\s*\\("

# Loop through each source file
for (file in source_files) {
  lines <- readLines(file)

  # Check if the file is an R script or a C++ source file
  if (grepl("\\.R$", file)) {
    # Search for R function declarations
    functions_in_file <- grep(r_function_pattern, lines, value = TRUE)

    # Clean function declarations to only include function names
    functions_in_file <- gsub("<-.*", "", functions_in_file)
    functions_in_file <- gsub("^.*\\$", "", functions_in_file)
  } else {
    # Search for C++ function declarations
    functions_in_file <- grep(cpp_function_pattern, lines, value = TRUE)

    # Clean function declarations to only include function names
    functions_in_file <- gsub("\\(.*", "", functions_in_file)
    functions_in_file <- gsub("^.*\\s+", "", functions_in_file)
  }

  function_names[[file]] <- functions_in_file
}

# Print the list of function names
function_names[sapply(function_names, length) != 0 ]
