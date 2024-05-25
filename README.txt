Solving Laplace Equations (David Baier):

Implemented the Laplace Beltrami Matrix, added further the vertex_weights_ to use the area for diagonal values. To compute the Laplace Belatrami matrix we iterate over all vertices and each neighbor per vertex. For each neighbor we get the vertex handle through the halfedge, the edge for the edge_weight. For each vertex the weight of the neighbors is summed. After all the neighbors were computed the summed value and the vertex_weights_ are added for the diagnoal value of the Laplace Beltrami matrix.
To get the missing values X we solve Ax = b with a Sparse Cholesky decomposistion SimplicialLDLT. 

The entries of X are the appended to the v_harmonic_function property by the vertex id.
There where some problems in verifying the function, since it was kind of not possible to change vertex0 and vertex1 manually (i guess do to some mismatch when clcking on the mesh.) Though the pattern between the default vertices seemed legit.
