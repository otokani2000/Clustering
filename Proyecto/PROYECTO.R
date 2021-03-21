#Se lee el archivo con la info. del blastp

tabla<-read.table("tabla_ORDENADA.txt")

tabla

#Tenemos 9992 alineamientos de 10000 que en teoría deberían ser. 
#Para buscar en que posiciones están los 8 faltantes de los 10000 hacemos lo sig.
#Guardamos la columna 1 y 2 en 2 vectores diferentes y vamos a ver cuales comparaciones o alineamientos no tienen valores de bistscores.
#1's para las casillas con valores, 0's en las casillas donde no se hayan producido alineamientos o comparaciones entre n y m secuencia. 

pos1<-factor(tabla[1][[1]])
pos2<-factor(tabla[2][[1]])

#Verificamos que longitud sea 10000
length(table(pos1,pos2)==0)

#Guardamos los índices de  las secuencias n y m donde no se produjeron alineamientos.
index<-which(table(pos1,pos2)==0)
index
class(which(table(pos1,pos2)==0))

#Guardamos los bitscores en un vector.

bitscores<-tabla[[3]]
length(bitscores)

#Al normalizar los valores de nuestra matriz va a tener valores entre 0-1.
#El valor 1 va a ser para los bitscores más altos y para valores en la diagonal.
#En los indices donde no se produjo alineamientos solo vamos a agregar 0's para poder generar posteriormente un matriz cuadrada.

contador<-0
for (x in index){
  bitscores<-append(bitscores,0,after=(x+contador-1))
  contador<-contador+1
}

#Verificamos que sean 10000 valores
length(bitscores)

#El valor de bitscore más alto se ecnuentra en la pos. 6162 y es 324
bitscores[6162]

#Gracias a que tenemos los alineamientos ordenados alfabeticamente por los nombres de las secuencias
#y sus respetivos bistscores ordenados, generamos una matriz cuadrada rellenando por filas. 
MATRIZ<-matrix(data=bitscores,ncol=100,nrow=100,byrow = TRUE)

# Agregamos los nombres para las filas y verificamos las dimensiones.
rownames(MATRIZ)<-levels(tabla[[2]][1])[1:100]
dim(MATRIZ)

#4.Normalizar las disimilitudes (d) para que queden en el rango [0,1].

#Se divide toda la matriz entre el bitscore más alto.
MATRIZ<-MATRIZ/324

#Rellenar 1's en la diagonal.
for (i in 1:100){
  MATRIZ[i,i]<-1  
  
}
#Tenemos una matriz de simiitud y para generar una de disimilitud hacemos lo sig:
MATRIZ<-1-MATRIZ

#El valor máximo debe ser 1 y el mínimo debe ser 0.
max(MATRIZ)
min(MATRIZ)
MATRIZ

#CLUSTERING

library(cluster)
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(ape))

#Leemos los datos

InputData  <- MATRIZ


#5. Correr clustering jerárquico y correr varios métodos para obtener el número de clusters.
#Run the hierarchical clustering and plot the dendogram

ccom <- hclust(dist(InputData), method = "ward.D2")
plot (ccom, hang = -1)

my_tree <- as.phylo(ccom)
write.tree(phy=my_tree, file="jerarquico.tree")

#Fijamos un área para graficar las 4 gráficas al mismo tiempo.

par(mfrow = c(1, 4))


#Fijamos los margenes para graficar

par(mar = c(2, 2, 2, 1) + 0.1)


#Construimos los dendogramas de diferentes métodos.

csin <- hclust(dist(InputData, method = "euclidean"), method = "single")
cave <- hclust(dist(InputData, method = "euclidean"), method = "average")
ccom <- hclust(dist(InputData, method = "euclidean"), method = "complete")
cwar <- hclust(dist(InputData, method = "euclidean"), method = "ward.D2")

#Generamos los coeficientes de clustering para cada método.
coeff_csin<-coef(csin)
coeff_cave<-coef(cave)
coeff_ccom<-coef(ccom)
coeff_cwar<-coef(cwar)

coeff_csin
coeff_cave
coeff_ccom
coeff_cwar

#6. Salvar el dendograma como árbol filogenético en formato Newick en R.

#Salvamos todos los trees de los diferentes métodos en formato Newick.
my_tree <- as.phylo(csin)
write.tree(phy=my_tree, file="csin.tree")

my_tree <- as.phylo(cave)
write.tree(phy=my_tree, file="cave.tree")

my_tree <- as.phylo(ccom)
write.tree(phy=my_tree, file="ccom.tree")

my_tree <- as.phylo(cwar)
write.tree(phy=my_tree, file="cwar.tree")

# Graficamos todos los dendogramas.

plot (csin, hang = -1, main = "Single")
rect.hclust(csin, k=20,  border=1:16)
csin20 <- cutree(csin, k=20)

plot (cave, hang = -1, main = "Average")
rect.hclust(cave, k=20,  border=1:16)
cave20 <- cutree(cave, k=20)

plot (ccom, hang = -1, main = "Complete")
rect.hclust(ccom, k=20,  border=1:16)
ccom20 <- cutree(ccom, k=20)

plot (cwar, hang = -1, main = "Ward.D")
rect.hclust(cwar, k=20,  border=1:16)
cwar20 <- cutree(cwar, k=20)


# Visualizamos el cluster ccom de otra manera 

cls3 <- cutree(ccom, k=4)
plot(InputData, xlim=c(0,8), ylim=c(0,8), col=cls3)
fviz_cluster(list(data = InputData, cluster = cls3))


#Ahora visualizamos el cluster csin.

csin <- hclust(dist(InputData), method = "single")
plot (csin, hang = -1)
rect.hclust(csin, k=4, border=2:4)

dend1 <- as.dendrogram (ccom)
dend2 <- as.dendrogram (csin)

dend_list <- dendlist(dend1, dend2)
tanglegram(dend1, dend2, main = paste("Entanglement =", round(entanglement(dend_list))))



# Determinamos el numero optimo de clusters.

#Methods: Total Within Sum of Squares (wss), silhouette, gap_stat
fviz_nbclust(InputData, FUN = hcut, method = "silhouette", k.max = 10, print.summary = TRUE)

#The silhouette method applied to hclust
fviz_nbclust(InputData, FUN = hcut, hc_func = "hclust", hc_method = "ward.D2", method = "silhouette", k.max = 20) +
  labs(subtitle = "Silhouette method")

#The silhouette method applied to agnes
fviz_nbclust(InputData, FUN = hcut, hc_func = "agnes", hc_method = "ward.D2", method = "silhouette", k.max = 20) +
  labs(subtitle = "Silhouette method")

#The silhouette method applied to diana
fviz_nbclust(InputData, FUN = hcut, hc_func = "diana", hc_method = "ward.D2", method = "silhouette", k.max = 20) +
  labs(subtitle = "Silhouette method")


