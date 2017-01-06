# MStates

Contains experimental code for evaluating microstates hypothesis.

##Learning techniques

###1) Silhouette criterion was calculated for predefined number of microstates:
![Informed](https://github.com/VlastaKoudelka/MStates/blob/master/Results/Informed_algorithms.jpeg)

No local optima can be observed.

###2) Unsupervised clustering techniques:

![Dbscan](https://github.com/VlastaKoudelka/MStates/blob/master/Results/dbscan_no_clst.jpeg)

DBSCAN estimates a number of clusters based on two parameters: maximal distance in one neighborhood (here we normalized distance by mean distance between all samples) and minimum number of neighbours defining a core member.

![Meanshift](https://github.com/VlastaKoudelka/MStates/blob/master/Results/MeanShift_no_clst.jpeg)

Meanshift algorithm is mainly sensitive to its **bandwidth** hyperparameter, which is a width of typically Gaussian kernel used for mean estimation.

###3) Stochastic Neighbour Embedding

Low perplexity           |  High perplexity
:-------------------------:|:-------------------------:
![lowperplexity](https://github.com/VlastaKoudelka/MStates/blob/master/Results/t-SNE%20GFP_peaks_perp_5.jpeg)   |  ![highperplexity](https://github.com/VlastaKoudelka/MStates/blob/master/Results/t-SNE%20GFP_peaks_perplexity30.jpeg)



The t-SNE maps all data samples from Ns=19 dimmensions to 2D space. It is particulary sensitive to perplexity parameter, which has the meaning of nuber of effective neighbours.

##PCA, Varimax, and P. M. methods (explained variance):
53% | 62% | 70%
:-------------------------:|:-------------------------:|:-------------------------:
![Pca](https://github.com/VlastaKoudelka/MStates/blob/master/Results/pca.jpg) | ![Varimax](https://github.com/VlastaKoudelka/MStates/blob/master/Results/varimax.jpg) | ![Pascual](https://github.com/VlastaKoudelka/MStates/blob/master/Results/pascual.jpg)


