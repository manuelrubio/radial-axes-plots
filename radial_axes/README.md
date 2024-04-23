# radial_axes

Code in Matlab, R, and Python for computing low-dimensional point representations of high-dimensional numerical data
according to several data visualization techniques based on a set of radial axes. The implemented methods are:

+ *'SC'*: **Star Coordinates**. A basic linear mapping.

       [1] KANDOGAN E.: Star coordinates: A multi-dimensional
         visualization technique with uniform treatment of dimensions. In
         Proceedings of the IEEE Information Visualization Symposium, 
         Late Breaking Hot Topics (2000), pp. 9-12.
		 
&nbsp;

+ *'RADVIZ'*: **RadViz**. Similar to Star Coordinates, but the (non-negative) data 
               is previously normalized so that sum of the rows of X is
               equal to 1. This generates a nonlinear mapping.

       [2] HOFFMAN P., GRINSTEIN G., MARX K., GROSSE I., STANLEY E.: DNA 
         visual and analytic data mining. In Proceedings of the 8th 
         conference on Visualization'97 (Los Alamitos, CA, USA, 1997), 
         VIS'97, IEEE Computer Society Press, pp. 437-441. 
        [doi:10.1109/VISUAL.1997.663916](https://ieeexplore.ieee.org/document/663916). 

       [3] RUBIO-SÁNCHEZ M., RAYA L., DÍAZ F., SANCHEZ A.: A comparative
         study between radviz and star coordinates. IEEE Transactions on 
         Visualization and Computer Graphics 22, 1 (Jan 2016), 619-628.
         [doi:10.1109/TVCG.2015.2467324](https://ieeexplore.ieee.org/document/7192699).
		 
&nbsp;

+ *'BIPLOT'*: **Principal component biplots**. A generalization of Principal
     Component Analysis.

       [4] GABRIEL K. R.: The biplot graphic display of matrices with 
         application to principal component analysis. Biometrika 58, 3 
         (Dec 1971), pp. 453--467. 
		 [doi:10.1093/biomet/58.3.453](https://academic.oup.com/biomet/article-abstract/58/3/453/233361).

&nbsp;

+ *'ARA'*: **Adaptable radial axes plots**. A hybrid approach between Star coordinates and biplots.

       [5] RUBIO-SÁNCHEZ M., SANCHEZ A., LEHMANN D. J.: Adaptable radial
         axes plots for improved multivariate data visualization.
         Computer Graphics Forum 36, 3 (2017), 389–399. 
         [doi:10.1111/cgf.13196](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.13196).

