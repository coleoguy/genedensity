

while (TRUE) {
  writeLines("This file was written to keep the SSD awake.", "../../keepawake.txt")
  file.remove("../../keepawake.txt")
  Sys.sleep(20)
}



