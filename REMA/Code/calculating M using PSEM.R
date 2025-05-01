# Calculating M using a phylogenetically-informed approach
# Code from https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Ffaf.12800&file=faf12800-sup-0001-DataS1.pdf

# To demonstrate the potential role for PCM in fisheries science, we re-analyze a foundational dataset compiled
# by Then et al.. We specifically download the file “Mlifehist_ver1.0.csv” and then include a copy as a data
# object in package phylosem to simplify the following demonstration. 

# Load packages
library(phylosem)
library(fishtree)
library(phylobase)
library(sem)

# Download tree
out = fishtree_complete_phylogeny()
tree = out[[1]]

# Load data object
data( Mlifehist_ver1_0 )
Data = Mlifehist_ver1_0

Data.in <- Mlifehist_ver1_0
filter(Data.in, Genus == "Sebastes" & Species == "ruberrimus")

# Reformat to match tree$tip.label
Data$Genus_species = factor( paste0(Data$Genus, "_", Data$Species) )

# Drop duplicates ... not dealing with variation among stocks within species
Data = Data[match(unique(Data$Genus_species),Data$Genus_species), ]

# log-transform to simplify later syuntax
Data = cbind( Data, "logM" = log(Data[,'M']),
              "logK" = log(Data[,'K']),
              "logtmax" = log(Data[,'tmax']),
              "logLinf" = log(Data[,'Linf']) )

# Identify species in both datasets
species_to_use = intersect( tree$tip.label, Data$Genus_species )

species_to_drop = setdiff( Data$Genus_species, tree$tip.label )

# Drop tips not present in trait-data
# Not strictly necessary, but helpful to simplify later plots
tree = ape::keep.tip( tree, tip=species_to_use )

# Drop trait-data not in phylogeny
# Necessary to define correlation among data
rows_to_use = which( Data$Genus_species %in% species_to_use )
Data = Data[rows_to_use,]

# Only include modeled variables in trait-data passed to phylosem
rownames(Data) = Data$Genus_species
Data = Data[,c('logM','logK','logtmax','logLinf')]


# Specify SEM structure
sem_structure = "
logK -> logtmax, b1
logLinf -> logtmax, b2
logtmax -> logM, a
"
# Grid-search model selection using AIC for transformations
Grid = expand.grid( "OU" = c(FALSE,TRUE),
                    "lambda" = c(FALSE,TRUE),
                    "kappa" = c(FALSE,TRUE) )
psem_grid = NULL
for( i in 1:nrow(Grid)){psem_grid[[i]] = phylosem( data=Data, tree = tree, sem = sem_structure, estimate_ou = Grid[i,'OU'], estimate_lambda = Grid[i,'lambda'], estimate_kappa = Grid[i,'kappa'], quiet = TRUE )}

for( i in 1:nrow(Grid)){psem_grid[[i]] = phylosem( data=Data, tree = tree, sem = sem_structure, estimate_ou = Grid[i,'OU'], estimate_lambda = Grid[i,'lambda'], estimate_kappa = Grid[i,'kappa'])}

# Extract AIC for each model and rank-order by parsimony
Grid$AIC = sapply( psem_grid, \(m) m$opt$AIC )
Grid = Grid[order(Grid$AIC,decreasing=FALSE),]

# Select model with lowest AIC
psem_best = psem_grid[[as.numeric(rownames(Grid[1,]))]]

# Load for plotting
library(phylosignal)

# Plot path diagram
my_fitted_DAG = as_fitted_DAG(psem_best)
plot(my_fitted_DAG, type="color")

# Plot using phylobase
my_phylo4d = as_phylo4d( psem_best )
barplot(my_phylo4d)

# Total, direct, and indirect effects
my_sem = as_sem(psem_best)
effects(my_sem)

# how do we find model predictions for a single species?
options(max.print = 1000000)
my_phylo4d

# label node ancestor edge.length node.type         logM        logK     logtmax  logLinf
# Sebastes_ruberrimus  154      366   7.3279420       tip -4.051285073 -3.29954373  4.79579055 6.491634

exp(-4.051285073)
# this is just equal to the input value, 0.0174

my_phylo4d %>% filter(label == "Sebastes_ruberrimus")


summary(my_sem)

print(my_sem)

