#include "Parameterization.hh"

void Parameterization::solve_harmonic_function(const std::vector<int>& _constraint_indices,
        const std::string& _property_name) {
    // ------------- IMPLEMENT HERE ---------
    // Set up Laplace-Beltrami matrix of the mesh.
    // For the vertices for which the constraints are added,
    // replace the corresponding row of the system with the constraint
    // Store the result in vertex property v_harmonic_function.
    // ------------- IMPLEMENT HERE ---------
    calc_weights();
    int n = mesh_.n_vertices();

    Eigen::SparseMatrix<double> A(n,n);
    Eigen::MatrixXd b = Eigen::VectorXd::Zero(n);

    std::vector<Eigen::Triplet<double>> triplets;

    // ToDo Update Laplace Beltrami Matrix and Matrix b
    for(auto vertex : mesh_.vertices()){

        double tmp_result = 0.0;
        for(auto vertex_halfedge_iter = mesh_.voh_iter(vertex); vertex_halfedge_iter.is_valid(); vertex_halfedge_iter++){
            auto edge_handle = mesh_.edge_handle(*vertex_halfedge_iter);
            auto vertex_handle = mesh_.to_vertex_handle(*vertex_halfedge_iter);
            auto edge_weight =  mesh_.property(edge_weight_, edge_handle);
    
            tmp_result += edge_weight;
            // i not j case
            triplets.emplace_back(Eigen::Triplet<double>(vertex.idx(), vertex_handle.idx(), edge_weight));
        }
        auto vw = mesh_.property(vertex_weight_, vertex);
        // i=j case
        triplets.emplace_back(Eigen::Triplet<double>(vertex.idx(), vertex.idx(), 1/vw - tmp_result));
    }


    A.setFromTriplets(triplets.begin(), triplets.end());
    //Iterate over all constained vertices and set them accordingly in A and b.
    int j = 0;
    for(auto indices : _constraint_indices){
        int ind = indices;
        if(indices > n)
            ind = j*(n-1);

        for(int i = 0; i < n; i++)
            A.coeffRef(ind, i) = 0;
        A.coeffRef(ind,ind) = 1;
    
        b(ind) = 1.0;
        ++j;
    }
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
    Eigen::VectorXd x = solver.solve(b);
    if(solver.info() != Eigen::Success)
        std::cout << "Error during solving" << std::endl;

    // --------------------------------------
    OpenMesh::VPropHandleT<double> v_harmonic_function;
    if(!mesh_.get_property_handle(v_harmonic_function, _property_name)) {
        mesh_.add_property(v_harmonic_function, _property_name);
        mesh_.property(v_harmonic_function).set_persistent(true);
     }

    // ------------- STORE RESULT HERE ---------
    for(auto vh:mesh_.vertices()){
        mesh_.property(v_harmonic_function, vh) = x(vh.idx());
    }
    // --------------------------------------

    color_coding(v_harmonic_function);
}

std::pair<int, int> Parameterization::get_intervals_borders(const double _iso0, const double _iso1,
        double _lower_bound, double _interval_size) {

    std::pair<int, int> intervals_borders;

    // ------------- IMPLEMENT HERE ---------
    // Given the values of the harmonic function that correspond to the two vertices in a triangle,
    // find the first and the last interval border that fall between the isovalues at the two vertices
    // Use std::pair to return the indices of the first and the last interval border.


    // --------------------------------------

    return intervals_borders;
}

void Parameterization::add_isoline_segment(const std::pair<int, int>& _borders01, const std::pair<int, int>& _borders02,
        const double _iso0, const double _iso1, const double _iso2, const Point& _p0, const Point& _p1, const Point& _p2,
        double _lower_bound, double _interval_size) {

    // ------------- IMPLEMENT HERE ---------
    // For each two edges of a triangle check if they are intersected by the same isoline.
    // If this is the case, compute the intersections using linear interpolation of the isovalues.
    // Add an isoline segment when the isoline indices for the two edges coincide
    // (segment_points_.push_back(p0); segment_points_.push_back(p1);)


    // --------------------------------------
}

void Parameterization::compute_isolines(TetrahedralMesh& _line_mesh, const std::vector<int>& _constraint_indices,
        const std::string& _property_name, int _num_intervals) {
    segment_points_.clear();

    OpenMesh::VPropHandleT<double> v_harmonic_function;
    if(mesh_.get_property_handle(v_harmonic_function, _property_name)) {
        auto lower_bound = mesh_.property(v_harmonic_function, OpenMesh::VertexHandle(_constraint_indices[1]));
        auto upper_bound = mesh_.property(v_harmonic_function, OpenMesh::VertexHandle(_constraint_indices[0]));
        auto interval_size = (upper_bound - lower_bound) / (double)_num_intervals;

        for (auto fh : mesh_.faces()) {
            std::vector<Point> pts;
            std::vector<double> iso_vals;
            for (auto fv_it = mesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
                iso_vals.push_back(mesh_.property(v_harmonic_function, *fv_it));
                pts.push_back(mesh_.point(*fv_it));
            }

            std::pair<int, int> borders01 = get_intervals_borders(iso_vals[0], iso_vals[1], lower_bound, interval_size);
            std::pair<int, int> borders12 = get_intervals_borders(iso_vals[1], iso_vals[2], lower_bound, interval_size);
            std::pair<int, int> borders02 = get_intervals_borders(iso_vals[0], iso_vals[2], lower_bound, interval_size);

            add_isoline_segment(borders01, borders02, iso_vals[0], iso_vals[1], iso_vals[2], pts[0], pts[1], pts[2], lower_bound, interval_size);
            add_isoline_segment(borders01, borders12, iso_vals[1], iso_vals[0], iso_vals[2], pts[1], pts[0], pts[2], lower_bound, interval_size);
            add_isoline_segment(borders02, borders12, iso_vals[2], iso_vals[0], iso_vals[1], pts[2], pts[0], pts[1], lower_bound, interval_size);
        }
    }

    //create mesh
    std::vector<OpenVolumeMesh::VertexHandle> vhs;
    for(auto pt : segment_points_)
        vhs.push_back(_line_mesh.add_vertex(pt));

    for(size_t i=0; i<vhs.size()/2; ++i)
        _line_mesh.add_edge(vhs[2*i], vhs[2*i + 1]);
}
