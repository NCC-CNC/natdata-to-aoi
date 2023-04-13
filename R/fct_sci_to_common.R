ch_cosewicid_to_name <- function(LUT, cosewicid, name_type) {
  x <- LUT %>% filter(COSEWIC_ID == cosewicid)
  if (name_type == "common") {
    name <- x$CommName_E[1]
  } else {
    name <- x$SciName[1]
  }
}

sar_cosewicid_to_name <- function(LUT, cosewicid, name_type) {
  x <- LUT %>% filter(COSEWICID == cosewicid)
  if (name_type == "common") {
    name <- x$COM_NAME_E[1]
  } else {
    name <- x$SCI_NAME[1]
  }
}

nsc_end_to_name <- function(LUT, sci_name) {
  x <- LUT %>% filter(NSC_SCI_NAME == sci_name)
  name <- str_to_title(x$COMMON_NAME)
  if (identical(name, character(0))) {
    name <- sci_name
  }
  if (is.na(name)) {
    name <- sci_name
  }
  return(name)  
}

nsc_sar_to_name <- function(LUT, sci_name) {
  x <- LUT %>% filter(CAN_SCI_NAME == sci_name)
  name <- str_to_title(x$CAN_COM_NAME)
  if (identical(name, character(0))) {
    name <- sci_name
  }
  if (is.na(name)) {
    name <- sci_name
  }
  return(name)
}

nsc_spp_to_name <- function(LUT, sci_name) {
  x <- LUT %>% filter(NATIONAL_SCIENTIFIC_NAME == sci_name)
  if (nrow(x) > 1) {
    x <- x[1,]
  }
  name <- str_to_title(x$NATIONAL_ENGL_NAME)
  
  if (identical(name, character(0))) {
    name <- sci_name
  }   
  if (is.na(name)) {
    name <- str_to_title(x$ENGLISH_COSEWIC_COM_NAME)
    if(is.na(name)) {
      name <- sci_name
    }
  }
  return(name)
}

iucn_to_name <- function(LUT, file_name) {
  x <- LUT %>% filter(File_Name == file_name)
  if (is.na(x$Common_Name)) {
    name <- x$Sci_Name
    if (!is.na(x$Season)) {
      name <- paste(x$Sci_Name, x$Season)
    }
  } else {
    name <- x$Common_Name
    if (!is.na(x$Season)) {
      name <- paste(x$Common_Name, x$Season)
    }
  }
  return(name)
}
