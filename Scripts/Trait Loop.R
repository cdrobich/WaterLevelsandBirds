
BSp<-read.csv("Data/Bird_community_composition.csv", row.names=1)
BTrt<-read.csv("Data/Traits_betweenyears_28MAY17_nmrwbl.csv", row.names=1)


library(vegan)

BSp<-BSp[,4:26] 

rownames(BTrt)
colnames(BSp)


# loop variables
sites <- row.names(BSp)
birds <- colnames(BSp)
traits <- colnames(BTrt)

# output table
Sum <- data.frame(matrix(as.numeric(0),ncol = length(traits), nrow = length(sites)))
rownames(Sum) <- sites
colnames(Sum) <- traits


# for loop
for (s in 1:length(sites)){                                       # Iterate across all rows and
  for (b in 1:length(birds)){                                     # columns in AB and check if
    if(BSp[row.names(BSp)==sites[s],names(BSp)==birds[b]]>0){                 # there is a count > 0 for a particular site and bird. If so,
      for (t in 1:length(traits)){                                  # check all the traits for
        if (birds[b] %in% row.names(BTrt)){                           # whether that bird is listed in the traits table and
          if(BTrt[birds[b],traits[t]]==1){                            # if that bird has trait (t). If so,
            Sum[row.names(Sum)==sites[s],names(Sum)==traits[t]]<-   # populate the column in Sum corresponding to that site and that trait with the sum of
              Sum[row.names(Sum)==sites[s],names(Sum)==traits[t]]+  # the count of that bird in AB and the count that
              BSp[row.names(BSp)==sites[s],names(BSp)==birds[b]]}               # was there previously.
        } else {print(paste("Skipped",birds[b],sep=" "))}           # this prints out the birds that are listed in AB but not in TR. Output = 'Merlin' and 'GRHE'
      }}}}



# output

write.csv(Sum, file = "Data/traits_matrix_Routput.csv")


