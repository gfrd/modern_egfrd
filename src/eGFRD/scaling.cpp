//
//
#include "scaling.hpp"
//
//inline double find_maximal_cylinder_height_to_plane(Vector3 start1, Vector3 unit_z, double height_static,
//                                                    double scaling_factor, PlanarSurface* plane)
//{
//    auto start2 = plane->project_point(start1).first;
//    auto w0 = start2 - start1;
//    auto upper_limit = w0.length();
//
//    auto f = [start1, unit_z, height_static, scaling_factor, plane, start2](double height_dynamic) {
//        Vector3 end1 = start1 + unit_z * (height_static + height_dynamic);
//        auto end2 = plane->project_point(end1).first;
//        auto squared_distance = squared_distance::line_segments(start1, end1, start2, end2);
//
//        // needed_distance is the cylinder radius required at given height, as determined by the scaling factor
//        // that ensures isotropic diffusion of IV
//        auto needed_distance = height_dynamic * scaling_factor;
//        return squared_distance - needed_distance*needed_distance;
//    };
//
//    if(f(upper_limit) > 0.0) {
//        // Even maximal height does not reach the plane, so return it directly
//        return upper_limit;
//    }
//
//    gsl_lambda<decltype(f)> F(f);
//    root_fsolver_wrapper solver;
//    solver.set(&F, 0, upper_limit);
//    try {
//        return solver.findRoot(1e-18, 1e-12, "scaling::find_maximal_cylinder_height_to_plane");
//    }
//    catch(const gfrd_exception& e)
//    {
//        Log("GFRD").warn() << "Exception whilst finding root of plane-cylinder scaling distance function";
//        return 0.0;
//    }
//}
//
//inline double find_maximal_cylinder_height_to_segment(Vector3 start1, Vector3 unit_z, double height_static,
//                                                      double scaling_factor, Vector3 start2, Vector3 end2,
//                                                      double other_radius)
//{
//    auto w0 = start2 - start1;
//    auto upper_limit = w0.length();
//
//    auto f = [start1, unit_z, height_static, scaling_factor, start2, end2, other_radius](double height_dynamic) {
//        Vector3 end1 = start1 + unit_z * (height_static + height_dynamic);
//        auto squared_distance = 0.0;
//
//        if(start2 != end2)
//        {
//            squared_distance = squared_distance::line_segments(start1, end1, start2, end2);
//        } else
//        {
//            // Start1 == end1, so the second line segment is actually a point
//            squared_distance = squared_distance::line_segment_point(start1, end1, start2);
//        }
//
//        // Needed_radius is the cylinder radius required at given height, as determined by the scaling factor
//        // that ensures isotropic diffusion of IV
//        auto needed_distance = height_dynamic * scaling_factor + other_radius;
//        return squared_distance - needed_distance*needed_distance;
//    };
//
//    if(f(upper_limit) > 0.0) {
//        // Even maximal height does not reach the plane, so return it directly
//        return upper_limit;
//    }
//
//    gsl_lambda<decltype(f)> F(f);
//    root_fsolver_wrapper solver;
//    solver.set(&F, 0, upper_limit);
//    try {
//        return solver.findRoot(1e-18, 1e-12, "scaling::find_maximal_cylinder_height_to_segment");
//    }
//    catch(const gfrd_exception& e)
//    {
//        Log("GFRD").warn() << "Exception whilst finding root of segment-cylinder scaling distance function";
//        return 0.0;
//    }
//}
//
//template<typename TCollector>
//struct max_height_collector
//{
//    max_height_collector(Vector3 base_position, Vector3 unit_z, double height_static, double scaling_factor) :
//            base_position(base_position), unit_z(unit_z), height_static(height_static), scaling_factor(scaling_factor), max_height(0.0) {}
//
//    void operator()(typename TCollector::const_iterator i, const Vector3&)
//    {
//        auto shell = (*i).second;
//        Vector3 start2, end2;
//        double other_radius = 0.0;
//        if(shell.shape() == Shell::Shape::SPHERE)
//        {
//            auto sphere = shell.get_sphere();
//            start2 = sphere.position();
//            end2 = sphere.position();
//            other_radius = sphere.radius();
//        }
//        else if (shell.shape() == Shell::Shape::CYLINDER)
//        {
//            auto cylinder = shell.get_cylinder();
//            start2 = cylinder.position() - cylinder.half_length() * cylinder.unit_z();
//            end2 = cylinder.position() + cylinder.half_length() * cylinder.unit_z();
//            other_radius = cylinder.radius();
//        }
//
//        double height = find_maximal_cylinder_height_to_segment(base_position, unit_z, height_static, scaling_factor, start2, end2, other_radius);
//        max_height = std::min(height, max_height);
//    }
//
//    double max_height;
//    Vector3 base_position;
//    Vector3 unit_z;
//    double height_static;
//    double scaling_factor;
//};
//
///**
// * Finds a maximal dynamical cylinder height (where total cylinder height is static height + dynamical height)
// * for a given cylinder, such that it does not touch any neighbouring structures or domains.
// *
// * In case of a Mixed2D3DPair or SinglePlanarInteraction, the static height is the distance from the 3D-diffusing
// * particle to their PlanarSurface + the length by which the cylinder exceeds beyond the other side of the PlanarSurface,
// * and the dynamic height is the length from the 3D-diffusing particle to the cylinder cap away from the PlanarSurface.
// *
// * @tparam TMatrixSpace
// * @param base_pos Position of the bottom of the cylinder. That is, the side of the 2D-diffusing particle.
// * @param unit_z Unit vector pointing from the bottom to the top of the cylinder.
// * @param height_static The static length of the cylinder, as described above.
// * @param scaling_factor The factor such that height_dynamic = radius / scaling_factor, that ensures isotropic diffusion of the IV.
// * @param matrixSpace The matrix space that divides the total world region into cells.
// * @param structures All structures in the world region.
// * @param ignored_structures A list of structure IDs that will not be compared against (e.g. attached PlanarSurfaces, and the world).
// * @return The maximal dynamical height of the cylinder, such that no neighbouring structures or domains are touched.
// */
//template <typename TMatrixSpace>
//double scaling::find_maximal_cylinder_height(Vector3 base_pos, Vector3 unit_z, double height_static,
//                                             double scaling_factor, const TMatrixSpace& matrixSpace,
//                                             const StructureContainer::structures_range structures,
//                                             const std::vector<StructureID>& ignored_structures)
//{
//    double max_dynamic_height = 10e6;
//
//    for(const auto& structure : structures)
//    {
//        if (std::find(ignored_structures.begin(), ignored_structures.end(), structure->id()) == ignored_structures.end()) {
//            auto plane = dynamic_cast<PlanarSurface *>(structure.get());
//            if (plane == nullptr) {
//                THROW_EXCEPTION(not_implemented,
//                                "Scaling cylinders up to structures other than PlanarSurfaces is not yet implemented.");
//            }
//
//            // TODO: add cyclic transpose to vectors
//
//            auto calculated_max_height = find_maximal_cylinder_height_to_plane(base_pos, unit_z,
//                                                                               height_static,
//                                                                               scaling_factor, plane);
//            max_dynamic_height = std::min(max_dynamic_height, calculated_max_height);
//        }
//    }
//
//    max_height_collector<TMatrixSpace> col(base_pos, unit_z, height_static, scaling_factor);
//    matrixSpace.each_neighbor_cyclic(matrixSpace.index(base_pos), col);
//    max_dynamic_height = std::min(max_dynamic_height, col.max_height);
//
//    return max_dynamic_height + height_static;
//}
//
//template<typename TCollector>
//struct max_radius_collector
//{
//    max_radius_collector(Vector3 start1, Vector3 end1) :
//            start1(start1), end1(end1), max_radius(1e6) {}
//
//    void operator()(typename TCollector::const_iterator i, const Vector3&)
//    {
//        auto shell = (*i).second;
//        Vector3 start2, end2;
//        double other_radius = 0.0;
//        if(shell.shape() == Shell::Shape::SPHERE)
//        {
//            auto sphere = shell.get_sphere();
//            other_radius = sphere.radius();
//            auto distance = std::sqrt(squared_distance::line_segment_point(start1, end1, sphere.position()));
//
//            if (distance == 0.0)
//            {
//                return;
//            }
//
//            max_radius = std::min(max_radius, distance - other_radius);
//        }
//        else if (shell.shape() == Shell::Shape::CYLINDER)
//        {
//            auto cylinder = shell.get_cylinder();
//            start2 = cylinder.position() - cylinder.half_length() * cylinder.unit_z();
//            end2 = cylinder.position() + cylinder.half_length() * cylinder.unit_z();
//            other_radius = cylinder.radius();
//            auto distance = std::sqrt(squared_distance::line_segments(start1, end1, start2, end2));
//
//            if (distance == 0.0)
//            {
//                return;
//            }
//
//            max_radius = std::min(max_radius, distance - other_radius);
//        }
//    }
//    Vector3 start1;
//    Vector3 end1;
//    double max_radius;
//};
//
//template <typename TMatrixSpace>
//double scaling::find_maximal_cylinder_radius(Vector3 start1, Vector3 end1, const TMatrixSpace& matrixSpace,
//                                             const StructureContainer::structures_range structures,
//                                             const std::vector<StructureID>& ignored_structures) {
//    double max_radius = 10e6;
//
//    for (const auto &structure : structures) {
//        if (std::find(ignored_structures.begin(), ignored_structures.end(), structure->id()) ==
//            ignored_structures.end()) {
//            auto plane = dynamic_cast<PlanarSurface *>(structure.get());
//            if (plane == nullptr) {
//                THROW_EXCEPTION(not_implemented,
//                                "Scaling cylinders up to structures other than PlanarSurfaces is not yet implemented.");
//            }
//
//            // TODO: add cyclic transpose to vectors
//
//            auto start2 = plane->project_point(start1).first;
//            auto end2 = plane->project_point(end1).first;
//            auto distance = std::sqrt(squared_distance::line_segments(start1, end1, start2, end2));
//            max_radius = std::min(max_radius, distance);
//        }
//    }
//
//    max_radius_collector<TMatrixSpace> col(start1, end1);
//    matrixSpace.each_neighbor_cyclic(matrixSpace.index(start1), col);
//    max_radius = std::min(max_radius, col.max_radius);
//
//    return max_radius;
//}


/***
 * Finds the minimal distance from a point in a plane to one of its edges, either in its X or Y directions.
 * If the point is outside of the plane, a distance of 0.0 is returned.
 * @param pos Vector in the plane
 * @param plane_id StructureID of the plane
 * @param world The world object the structure resides in
 * @return double, the minimal orthogonal distance to an edge of the plane
 */
double scaling::dist_to_plane_edge(const Vector3 pos, const StructureID plane_id, const World& world)
{
    auto structure = world.get_structure(plane_id);
    auto plane = dynamic_cast<PlanarSurface*>(structure.get());
    if (plane == nullptr) {
        // Structure is not a plane
        THROW_EXCEPTION(illegal_argument, "Structure given in dist_to_plane_edge() is not a plane");
    }

    auto projected = plane->project_point(pos).first;
    auto half_extent = plane->shape().half_extent();
    auto plane_pos = plane->shape().position();
    auto offset_x = fabs(projected.X() - plane_pos.X());
    auto offset_y = fabs(projected.Y() - plane_pos.Z()); // Z-axis and X-axis form a basis for the plane
    auto distance = std::max(0.0, std::min(half_extent.X() - offset_x, half_extent.Y() - offset_y));
    return distance;
}