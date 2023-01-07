library(dplyr)

d <- readxl::read_xlsx("~/Downloads/sdmTMB 2023 TESA workshop participants.xlsx", sheet = 2)
names(d) <- tolower(names(d))
glimpse(d)

d[, 2] <- NULL
names(d)[7] <- "regional_priority"
names(d)[8] <- "location"
names(d)[3] <- "first"
names(d)[2] <- "last"
d$`supervisor/program` <- NULL

glimpse(d)
nrow(d)

d$email <- tolower(d$email)
sort(unique(d$email))

domains <- strsplit(d$email, "@")
domains <- unlist(lapply(domains, function(x) x[[2]]))
unique(domains)
table(domains)

d$email_formatted <- paste0(d$first, " ", d$last, " <", d$email, ">")
d <- arrange(d, last)
# View(d)

dups <- d[duplicated(select(d, last, first)), ]
# View(dups)

d <- d[!duplicated(select(d, last, first)), ]
nrow(d)

email_list <- paste(d$email_formatted, collapse = "; ")
print(email_list)

# already added:
already <- scan("~/Desktop/emailed.txt", what = "character")
already <- gsub(";", "", already)
already <- tolower(already)

to_add <- !d$email %in% already
paste(d$email[to_add], collapse = "; ")

already[!already %in% d$email]

nrow(d)

# 80
# 8-9
# 6
# 4-5?
# designate someone to count and declare how many groups and leaders
# I'll click the magic button

# helpers:
# - me,
# - philina, eric, lewis, lindsay, cam, semra, paul
#
# - patrick(?), jessica(?), ...
# - extras: Jillian, Dan, Catarina, Mackenzie, Noel, Elise, ...
# 8-10

#
