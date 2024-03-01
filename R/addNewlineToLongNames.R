#' Add Newline to Long Names in a Data Frame Column
#'
#' This function takes a dataframe and a column name, then modifies the contents of that column
#' by inserting newline characters into long variable names. The aim is to ensure that no line
#' exceeds a specified maximum character length, improving readability in plots or tables where
#' space is limited.
#'
#' @param data A dataframe containing the column to be modified.
#' @param column_name The name of the column within `data` that contains the variable names to be adjusted.
#' @param max_char The maximum number of characters allowed per line before inserting a newline.
#'        Defaults to 30 characters.
#'
#' @return The modified dataframe with newline characters added to the specified column's entries.
#'
#' @examples
#' df <- data.frame(variable_name = c("This is a very long variable name that needs to be split",
#'                                    "Short name",
#'                                    "Another extremely long variable name that requires splitting"))
#' df <- addNewlineToLongNames(df, "variable_name")
#' print(df)
addNewlineToLongNames <- function(data, column_name, max_char = 30) {
  # Apply the newline insertion logic to each element of the specified column
  data[[column_name]] <- sapply(data[[column_name]], function(name) {
    # Split the name into words
    words <- strsplit(name, " ")[[1]]
    new_name <- ""  # Initialize the new name
    current_line <- ""  # Keep track of the current line being constructed

    # Iterate through each word
    for (word in words) {
      # Check if adding the next word exceeds the maximum character limit
      if (nchar(current_line) + nchar(word) + nchar(ifelse(nchar(current_line) == 0, "", " ")) <= max_char) {
        # If not, add the word to the current line
        current_line <- paste0(current_line, ifelse(nchar(current_line) == 0, "", " "), word)
      } else {
        # If it does, add the current line to `new_name` and start a new line with the word
        new_name <- paste0(new_name, ifelse(nchar(new_name) == 0, "", "\n"), current_line)
        current_line <- word
      }
    }

    # Add the last line to `new_name`
    paste0(new_name, ifelse(nchar(new_name) == 0, "", "\n"), current_line)
  })

  return(data)
}
