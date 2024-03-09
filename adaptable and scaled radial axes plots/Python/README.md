# Calculate Mapping

&nbsp;

## Description

This function computes the mappings of several radial axes methods: Star Coordinates, RadViz and Adaptable Radial Axes Plots.

&nbsp;

## About the function

**mapping***(algorithm, X, V, W=None, vector_norm=None, chosen_variable=None)*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Returns an N by m matrix P with the low-dimensional embeddings stored.

#### Parameters

&nbsp;&nbsp;&nbsp;&nbsp; **algorithm: *string***. Defines the radial method. Posibles values: SC, RadViz, SRA, Adaptable, Adaptable exact or Adaptable ordered.


&nbsp;&nbsp;&nbsp;&nbsp; **X: *matrix***. Is an N by n matriz whose rows contains the n dimensional data samples.


&nbsp;&nbsp;&nbsp;&nbsp; **V: *matrix***. Is an n by m matrix whose rows defines the method's axis vectors.


&nbsp;&nbsp;&nbsp;&nbsp; **W: *matrix, optional***. Is an n by n diagonal matrix defining nonnegative weights for each variable.


&nbsp;&nbsp;&nbsp;&nbsp; **vector_norm: *Int or string, optional***. Is the vector norm associated with adaptable radial axes plots. Posibles values: 1, 2 or Inf.


&nbsp;&nbsp;&nbsp;&nbsp; **chosen_variable: *Int, optional***. Is the selected attribute for constrained adaptable radial axes plots.

#### Returns


&nbsp;&nbsp;&nbsp;&nbsp; **out: *matrix***. The low-dimensional embeddings are stored in the N by m matrix P.


&nbsp;

## Examples
	import numpy as np

	N = 100
	n = 5
	m = 2
	X = np.random.randn(N, n)
	V = np.random.randn(n, 2)
	W = np.identity(n)

	mapping('SC', X, V)
	mapping('RadViz', X, V)
	mapping('Adaptable', X, V, W, 1)
	mapping('Adaptable', X, V, W, 2)
	mapping('Adaptable', X, V, W, 'Inf')
	mapping('Adaptable exact', X, V, W, 1, 0)
	mapping('Adaptable exact', X, V, W, 2, 0)
	mapping('Adaptable exact', X, V, W, 'Inf', 0)
	mapping('Adaptable ordered', X, V, W, 1, 0)
	mapping('Adaptable ordered', X, V, W, 2, 0)
	mapping('Adaptable ordered', X, V, W, 'Inf', 0)

	
