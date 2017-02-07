# MStates

Contains experimental codes for evaluating microstates hypothesis.

##Learning techniques

Silhouette approach | DBSCAN | MeanShift
:--:|:--:|:--:
![Informed](https://github.com/VlastaKoudelka/MStates/blob/master/Results/Informed_algorithms.jpeg) | ![Dbscan](https://github.com/VlastaKoudelka/MStates/blob/master/Results/dbscan_no_clst.jpeg) | ![Meanshift](https://github.com/VlastaKoudelka/MStates/blob/master/Results/MeanShift_no_clst.jpeg)
Silhouette criterion was calculated for predefined number of microstates. No local optima can be observed. | DBSCAN estimates a number of clusters based on two parameters: **maximal distance** in one neighborhood (here we normalized distance by mean distance between all samples) and **minimum number of neighbours** defining a core member. | Meanshift algorithm is mainly sensitive to its **bandwidth** hyperparameter, which is a width of typically Gaussian kernel used for mean estimation (here we normalized distance by mean distance between all samples).

##Stochastic Neighbour Embedding

Low perplexity           |  High perplexity
:-------------------------:|:-------------------------:
![lowperplexity](https://github.com/VlastaKoudelka/MStates/blob/master/Results/t-SNE%20GFP_peaks_perp_5.jpeg)   |  ![highperplexity](https://github.com/VlastaKoudelka/MStates/blob/master/Results/t-SNE%20GFP_peaks_perplexity30.jpeg)



The t-SNE maps all data samples from Ns=19 dimmensions to 2D space. It is particulary sensitive to perplexity parameter, which has the meaning of nuber of effective neighbours.

##PCA, Varimax, and P. M. methods

###1) Spatial component analysis and P.M. algorithm comparison and explained variance
constr. 53% relaxed 88.5% |constr. 62% relaxed 88.5% | constr. 70% relaxed NA
:--:|:--:|:--:
![Pca](https://github.com/VlastaKoudelka/MStates/blob/master/Results/pca.jpg) | ![Varimax](https://github.com/VlastaKoudelka/MStates/blob/master/Results/varimax.jpg) | ![Pascual](https://github.com/VlastaKoudelka/MStates/blob/master/Results/pascual.jpg)
####PCA microstate spatial covariance matrix:

     1.0000    0.0000   -0.0000   -0.0000
     0.0000    1.0000    0.0000   -0.0000
    -0.0000    0.0000    1.0000   -0.0000
    -0.0000   -0.0000   -0.0000    1.0000

####Varimax microstate spatial covariance matrix:

     1.0000    0.2296   -0.4951   -0.4374
     0.2296    1.0000   -0.0190   -0.5923
    -0.4951   -0.0190    1.0000    0.4977
    -0.4374   -0.5923    0.4977    1.0000   

Notice that neither P.M. nor Varimax microstates are the orthogonal components. In this case, the microstates are projections of the time series to the components. That is the reason why Varimax microstates are not nescessary orthogonal (they are more similar to P.M. MSs from this point of view). 

####P.M. microstate spatial covariance matrix:
    
    1.0000    0.6895    0.5412    0.0263
    0.6895    1.0000    0.7157   -0.5854
    0.5412    0.7157    1.0000   -0.0007
    0.0263   -0.5854   -0.0007    1.0000
