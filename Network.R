library(tidyverse)
library(microViz)

# load in the following function
'%!in%' <- function(x,y) {
  !('%in%'(x,y))
}
# END
physeq_Lake22 = physeq_Lake22 %>% tax_fix() %>% tax_transform(trans = "hellinger", rank = "Genus")

OTU1 = as(otu_table(physeq_Lake22),"matrix")
as.data.frame(OTU1)
colnames(OTU1) <- as.data.frame(tax_table(physeq_Lake22))$Genus
as.matrix(OTU1)


spieceasi.net <- multi.spiec.easi(list(OTU1), method = 'mb',lambda.min.ratio = 1e-3,nlambda = 40,
                                  icov.select.params = list(rep.num = 50, thresh = 0.05, ncores=16))

#Create a simple plot using igraph
ig.mb     <- adj2igraph(getRefit(spieceasi.net))
am.coord <- layout.fruchterman.reingold(ig.mb)
plot(ig.mb, layout=am.coord, vertex.size=1, vertex.label=NA, main="MB")

#Extract the adjacency matrix from the spiec.easi object. This
#indicates which pairs of OTUs are adjacent or not in the graph

spieceasi.matrix <- symBeta(getOptBeta(spieceasi.net), mode='maxabs')
spieceasi.matrix.dsc <- spieceasi.matrix
spieceasi.matrix <- as.matrix(spieceasi.matrix)

#Add the OTU names to the adjacency matrix.
rownames(spieceasi.matrix) <- colnames(OTU1)
colnames(spieceasi.matrix) <- colnames(OTU1)
otu.names <- colnames(OTU1)
#Build a weighted network from the adjacency matrix. The
#edges in a weighted network represent the strength of association between OTUs.

net <- graph.adjacency(spieceasi.matrix, mode = "undirected",
                       weighted = TRUE, diag = FALSE)
V(net)$name <- otu.names                       

#Convert the edge weights into distances, where larger weights
#become shorter distances, and then output a distance-based network

net.dist <- net
max(abs(E(net.dist)$weight))
weights.dist <- 1 - abs(E(net.dist)$weight)
E(net.dist)$weight <- weights.dist                           

#Convert the weighted network to a separate absolute network
net.abs <- net
E(net.abs)$weight <- abs(E(net.abs)$weight) 


#Calculate centrality metrics and create a summary

# alpha centrality
net.alpha <- alpha.centrality(net)
# degree distribution
net.strength <- strength(net.abs)
# betweenness centrality
bet <- betweenness(net.dist,v = V(net.dist))
# make a summary of centrality metrics
summary_cent <- as.data.frame(net.alpha)
colnames(summary_cent) <- ("Alpha_centrality")
rownames(summary_cent) <- colnames(OTU1)
summary_cent$Weighted_vertex_degree <- net.strength
summary_cent$Betweenness_centrality <- bet
metrics <- summary_cent

wt <- cluster_louvain(net, weights = E(net.dist)$weight)
temp <- V(net)$name
temp <- as.data.frame(temp)
temp$louvain <- membership(wt)
V(net)$louvain <- temp$louvain                         

#See which nodes have been put into modules with less than or equal to three members and then combine them into a single
#group that should be consider as not having been assigned a module

length(unique(temp$louvain))
summary_modules <- data.frame(table(temp$louvain))
colnames(summary_modules) <- c("louvain", "n")
summary_modules
modules <- as.numeric(summary_modules$louvain[which(summary_modules$n
                                                    >3)])

x <- max(modules)+1
for (i in c(1:length(temp$temp))) {
  if(temp$louvain[i] %!in% modules){
    temp$louvain[i] <- paste(x)
  }}

modules <- temp

modules$louvain <- as.numeric(modules$louvain)

modules <- modules[order(modules$louvain),]

module.lookup <-
  data.frame("louvain"=unique(modules$louvain),"new_louvain" =
               c(1:length(unique(modules$louvain))))

new <- merge(modules,module.lookup)

modules <- new

modules <- modules[,2:3]

summary_modules <- data.frame(table(modules$new_louvain))

summary_modules

max(modules$new_louvain)



## export the network to Gephi to view the graph

# melt the network to prepare for gephi
spieceasi.matrix.m <- melt(spieceasi.matrix)
#name cols
colnames(spieceasi.matrix.m) <- c("source","target","weight")
# get the names of all nodes
node.names <-
  unique(c(as.character(unique(spieceasi.matrix.m$source)),as.character
           (unique(spieceasi.matrix.m$target))))
# number them as an alphabetical node list, write to a csv
node.names <- as.data.frame(node.names)
node.names <- as.data.frame(node.names)
node.names$node_number <- c(1:length(node.names$node.names))
node.names$node_number2 <- c(1:length(node.names$node.names))
colnames(node.names) <- c("Taxonomy", "Label", "Id")
row.names(node.names) <- node.names$Taxonomy
row.names(modules) <- modules$temp
modules <- modules[order(modules$temp),]
row.names(node.names) ==row.names(metrics)
row.names(node.names) ==row.names(modules)
node.names.final <- cbind(node.names, metrics,modules)
write.table(node.names.final, "node.names.csv", sep = ",", row.names
            = FALSE)


# convert node names to numbers, write to a csv
temp <- merge(x = spieceasi.matrix.m, y = node.names, by.x =
                "source", by.y = "Taxonomy")
# create the edge list
colnames(temp) <-
  c("source","target","weight","remove","source_number")
temp <- temp[,-4]
edge.list <- merge(x = temp, y = node.names, by.x = "target", by.y =
                     "Taxonomy")
colnames(edge.list) <- c("source","target","weight","source.number",
                         "target.number")
edge.list <- edge.list[,c(3,4,6)]
colnames(edge.list) <- c("weight","source","target")
edge.list$Type <- "Undirected"
negative <- ifelse(edge.list$weight<0, "negative", "positive")
edge.list$Negative <- negative
edge.list$weight <- abs(edge.list$weight)
edge.list <- edge.list[which(abs(edge.list$weight)>0),]
write.table(edge.list, "edge_list.csv", sep = ",", row.names = FALSE)





