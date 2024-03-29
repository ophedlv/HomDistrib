# Rscript ./plot_genes.R

# Import package
library("plyr")
library("dplyr")
library("stringr")
library("ggplot2")
library("reshape2")
library("maptools")
library("gpclib")


# clear all
rm(list = ls())

# Command line parameter : file
infile <- commandArgs(trailingOnly = TRUE)
path_root <- unlist(strsplit(infile, "res"))[1]

#---------------
# World map data
#---------------
shape <- readShapeSpatial(paste(path_root,"/data/Continents/continent.shp",sep="")) 
shape@data$id <- rownames(shape@data)
shape.fort <- fortify(shape, 'ID') 
shape.fort <- shape.fort[order(shape.fort$order), ] 

#--------------------------------------
# Abundance Matrix (geneID -> sampleID)
#--------------------------------------

# Path and file id
path_res <- dirname(infile)
name <- gsub("^(.*)sub_matrix_abund[.]txt$", "\\1", basename(infile))

# Create plot folder (if it does not exist)
plot_folder <- paste(path_res,"/plot",sep="")
if(!(dir.exists(plot_folder))) {dir.create(plot_folder)}


# Read filtered data
d <- read.table(infile,header=T, stringsAsFactor = F)

d_long <- melt(d, id.vars = c("geneID"), variable.name= "sampleID", value.name = "abund")
d_long$sampleID <- gsub("^.","",d_long$sampleID)
d_long <- filter(d_long, abund != 0)


#--------------------------------
# Correspondance (sampleID -> id)
#--------------------------------

# Read correspondance file
correspondance <- read.csv(paste(path_root,"/data/correspondance",sep=""), sep = "\t", header = T)

for(i in 1:dim(d_long)[1]) {
	flag = 0
	j = 1
	while(flag == 0) {
		if(correspondance$sampleID[j] == d_long$sampleID[i]) {
			d_long$stationID[i] = paste(correspondance$station[j],correspondance$depth[j],sep="")
			flag = 1
		}
		j = j+1
	}
}

#--------------------------------------------
# Environmental Data (stationid->coordinates)
#--------------------------------------------

# Read environmental data
env <- read.csv(paste(path_root,"/data/envdata20170215.tsv",sep=""), sep = "\t", stringsAsFactor = F) # Modify
env <- env[-seq(585,593),1:3] # remove unknown stations

# Formatting data
for (i in 1:length(env[,1])) {
	n <- nchar(env[i,1])
	env[i,1] <- substr(env[i,1], 2, n)
	env$prof[i] <- substr(env[i,1], n-3, n)
} 

# Coordinates extraction
for (i in 1:dim(d_long)[1]) {
	flag = 0
	j = 1
	while(flag == 0 & j < length(env[,1])) {
		if(d_long$stationID[i] == env[j,1]) {
			d_long$lat[i] = env$lat[j]
			d_long$lon[i] = env$lon[j]
			d_long$depth[i] = env$prof[j]
			flag = 1
		}
		j = j+1
	}
	if(flag == 0) {
		d_long$lat[i] = NA
		d_long$lon[i] = NA
		d_long$depth[i] = NA
	}
}

#-------------
# Plot
#-------------

# occ <- data.frame(table(d_long$geneID))
# occ <- occ[order(occ$Freq, decreasing = TRUE), ]
# genes <- occ$Var1[1:5] # n number of principal homologous
genes <- unique(d_long$geneID)

# Abundance of the n major homologous at each station
d_main <- data.frame()
for (i in 1:length(genes)) {
	tmp <- filter(d_long,geneID==genes[i])
	d_main <- rbind(d_main,tmp)
}
p <- ggplot() +
	geom_bar(aes(x = factor(geneID), fill = factor(depth)),data=d_main) +
	xlab("Gene ID") +
	ylab("Number station") +
	theme(text = element_text(size=5), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(name,"hist.pdf",sep=""), plot = p, path = plot_folder)


# Number of different homologous present at each station
d_group <- group_by(d_long, stationID, lat, lon, depth)
d_group <- data.frame(summarise(d_group, nb_gene = length(unique(geneID)), sum_abund = sum(abund), mean_abund = mean(abund, na.rm = T), sd_abund = sd(abund, na.rm = T), max_abund = max(abund), min_abund = min(abund)))
depth_value <- na.omit(unique(d_group$depth))
d_group$nb_gene <- as.factor(d_group$nb_gene)


for (i in 1:length(depth_value)) {
	p <- ggplot() +
		geom_polygon(aes(long, lat, group=group), data = shape.fort, colour='black',fill='white') + theme_bw() +
		geom_point(aes(x=longitude,y=latitude),data=env, colour = "grey") +
		geom_point(aes(x=lon, y=lat, colour = mean_abund, size = nb_gene),data=filter(d_group, depth == depth_value[i])) +
		scale_colour_gradient(low="yellow", high="red") +
		ggtitle(paste("Abondance moyenne",depth_value[i])) +
		coord_quickmap()
	ggsave(paste(name,paste(depth_value[i],"_map_mean.pdf",sep=""),sep=""), plot = p, path = plot_folder)

	p <- ggplot() +
		geom_polygon(aes(long, lat, group=group), data = shape.fort, colour='black',fill='white') + theme_bw() +
		geom_point(aes(x=longitude,y=latitude),data=env, colour = "grey") +
		geom_point(aes(x=lon, y=lat, colour = sum_abund, size = nb_gene),data=filter(d_group, depth == depth_value[i])) +
		scale_colour_gradient(low="yellow", high="red") +
		ggtitle(paste("Somme totale abondance",depth_value[i])) +
		coord_quickmap()
	ggsave(paste(name,paste(depth_value[i],"_map_sum.pdf",sep=""),sep=""), plot = p, path = plot_folder)

	p <- ggplot() +
		geom_polygon(aes(long, lat, group=group), data = shape.fort, colour='black',fill='white') + theme_bw() +
		geom_point(aes(x=longitude,y=latitude),data=env, colour = "grey") +
		geom_point(aes(x=lon, y=lat, colour = max_abund, size = nb_gene),data=filter(d_group, depth == depth_value[i])) +
		scale_colour_gradient(low="yellow", high="red") +
		ggtitle(paste("Abondance maximale",depth_value[i])) +
		coord_quickmap()
	ggsave(paste(name,paste(depth_value[i],"_map_max.pdf",sep=""),sep=""), plot = p, path = plot_folder)
		
	p <- ggplot() +
		geom_polygon(aes(long, lat, group=group), data = shape.fort, colour='black',fill='white') + theme_bw() +
		geom_point(aes(x=longitude,y=latitude),data=env, colour = "grey") +
		geom_point(aes(x=lon, y=lat, colour = min_abund, size = nb_gene),data=filter(d_group, depth == depth_value[i])) +
		scale_colour_gradient(low="yellow", high="red") +
		ggtitle(paste("Abondance minimale",depth_value[i])) +
		coord_quickmap()
	ggsave(paste(name,paste(depth_value[i],"_map_min.pdf",sep=""),sep=""), plot = p, path = plot_folder)
}

#------------
# Saving file
#------------

file_coord <- paste(name,"coord.csv",sep="")
write.csv(d_long,paste(path_res,file_coord,sep="/"),row.names=F)
file_group <- paste(name,"group.csv",sep="")
write.csv(d_group,paste(path_res,file_group,sep="/"),row.names=F)

