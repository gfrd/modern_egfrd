#pragma once

#include "Vector3.hpp"
#include "StructureContainer.hpp"
#include "Shell.hpp"
#include <helperFunctionsGf.hpp>
#include <Logger.hpp>
#include "helperFunctions.hpp"
#include "PlanarSurface.hpp"
#include "World.hpp"

namespace scaling
{
    double test();

    template <typename TMatrixSpace>
    double find_maximal_cylinder_height(Vector3 base_pos, Vector3 unit_z, double height_static,
                                                 double scaling_factor, const TMatrixSpace& matrixSpace,
                                                 const World& world,
                                                 const std::vector<StructureID>& ignored_structures,
                                                 const std::vector<ShellID>& ignored_shells,
                                                 bool social_scaling=true);

    template <typename TMatrixSpace>
    double find_maximal_cylinder_radius(Vector3 start1, Vector3 end1, const TMatrixSpace& matrixSpace,
                                                 const World& world,
                                                 const std::vector<StructureID>& ignored_structures,
                                                 const std::vector<ShellID>& ignored_shells,
                                                 bool social_scaling=true);
}



// -----------------



inline double social_correction(double size, const Shell::Code shell_type)
{
    switch(shell_type)
    {
        case Shell::Code::INIT:
            // Share available space equally with neighbouring init shell
            return size / 2;
        case Shell::Code::MULTI:
            // Give multi shell some extra space to move around
            return size / 1.1;
        case Shell::Code::NORMAL:
            // TODO: add a globally tweakable distance divider, that gets optimised in equilibration phase of simulation
            return size / 1.1;
        default:
            THROW_EXCEPTION(illegal_argument, "Shell is not of init, multi, or normal type");
            break;
    }
}


inline double find_maximal_cylinder_height_to_plane(Vector3 start1, Vector3 unit_z, double height_static,
                                                    double scaling_factor, PlanarSurface* plane)
{
    auto start2 = plane->project_point(start1).first;
    auto w0 = start2 - start1;
    auto upper_limit = w0.length();

    auto f = [start1, unit_z, height_static, scaling_factor, plane, start2](double height_dynamic) {
        Vector3 end1 = start1 + unit_z * (height_static + height_dynamic);
        auto end2 = plane->project_point(end1).first;
        auto squared_distance = squared_distance::line_segments(start1, end1, start2, end2);

        // needed_distance is the cylinder radius required at given height, as determined by the scaling factor
        // that ensures isotropic diffusion of IV
        auto needed_distance = height_dynamic * scaling_factor;
        return squared_distance - needed_distance*needed_distance;
    };

    if(f(upper_limit) > 0.0) {
        // Even maximal height does not reach the plane, so return it directly
        return upper_limit;
    }

    gsl_lambda<decltype(f)> F(f);
    root_fsolver_wrapper solver;
    solver.set(&F, 0, upper_limit);
    try {
        return solver.findRoot(1e-18, 1e-12, "scaling::find_maximal_cylinder_height_to_plane");
    }
    catch(const gfrd_exception& e)
    {
        Log("GFRD").warn() << "Exception whilst finding root of plane-cylinder scaling distance function";
        return 0.0;
    }
}

inline double find_maximal_cylinder_height_to_segment(Vector3 start1, Vector3 unit_z, double height_static,
                                                      double scaling_factor, Vector3 start2, Vector3 end2,
                                                      double other_radius)
{
    auto w0 = start2 - start1;
    auto upper_limit = w0.length();

    auto f = [start1, unit_z, height_static, scaling_factor, start2, end2, other_radius](double height_dynamic) {
        Vector3 end1 = start1 + unit_z * (height_static + height_dynamic);
        auto squared_distance = 0.0;

        if(start2 != end2)
        {
            squared_distance = squared_distance::line_segments(start1, end1, start2, end2);
        } else
        {
            // Start2 == end2, so the second line segment is actually a point
            squared_distance = squared_distance::line_segment_point(start1, end1, start2);
        }

        // Needed_distance is the cylinder radius required at given height, as determined by the scaling factor
        // that ensures isotropic diffusion of IV, plus the radius of the other object.
        auto needed_distance = height_dynamic * scaling_factor + other_radius;

        // When the squared distance and needed distance^2 are equal, this function passes through zero
        return squared_distance - needed_distance*needed_distance;
    };

    if(f(upper_limit) > 0.0) {
        // Even maximal height does not reach the plane, so return it directly
        return upper_limit;
    }

    gsl_lambda<decltype(f)> F(f);
    root_fsolver_wrapper solver;
    solver.set(&F, 0, upper_limit);
    try {
        return solver.findRoot(1e-18, 1e-12, "scaling::find_maximal_cylinder_height_to_segment");
    }
    catch(const gfrd_exception& e)
    {
        Log("GFRD").warn() << "Exception whilst finding root of segment-cylinder scaling distance function";
        return 0.0;
    }
}

template<typename TCollector>
struct max_height_collector
{
    max_height_collector(const World& world, Vector3 base_position, Vector3 unit_z, double height_static, double scaling_factor,
                         const std::vector<ShellID>& ignored_shells) :
            world(world), base_position(base_position), unit_z(unit_z), height_static(height_static),
            scaling_factor(scaling_factor), ignored_shells(ignored_shells), is_limited(false), max_height(1e6) {}

    void operator()(typename TCollector::const_iterator i, const Vector3&)
    {
        auto shellID = (*i).first;
        auto shell = (*i).second;

        if(std::find(ignored_shells.begin(),ignored_shells.end(), shellID) != ignored_shells.end())
        {
            // Skip shells in ignored list
            return;
        }

        Vector3 start2, end2;
        double other_radius = 0.0;
        if(shell.shape() == Shell::Shape::SPHERE)
        {
            auto sphere = shell.get_sphere();
            start2 = sphere.position();
            end2 = sphere.position();
            other_radius = sphere.radius();

            // Place the base closest to the sphere, given a cyclic world.
            base_position = world.cyclic_transpose(base_position, sphere.position());
        }
        else if (shell.shape() == Shell::Shape::CYLINDER)
        {
            auto cylinder = shell.get_cylinder();
            start2 = cylinder.position() - cylinder.half_length() * cylinder.unit_z();
            end2 = cylinder.position() + cylinder.half_length() * cylinder.unit_z();
            other_radius = cylinder.radius();

            // Place the base closest to the cylinder, given a cyclic world.
            base_position = world.cyclic_transpose(base_position, cylinder.position());
        }

        double height = find_maximal_cylinder_height_to_segment(base_position, unit_z, height_static, scaling_factor, start2, end2, other_radius);

        if(height < max_height)
        {
            max_height = height;
            is_limited = true;
            limiting_shell_type = shell.code();
        }
    }

    const World& world;
    double max_height;
    Vector3 base_position;
    Vector3 unit_z;
    double height_static;
    double scaling_factor;
    const std::vector<ShellID>& ignored_shells;
    bool is_limited; // true if the height is limited by a neighbouring shell
    Shell::Code limiting_shell_type; // shell type of neighbouring cell that limits height
};

template<typename TCollector>
struct max_radius_collector
{
    max_radius_collector(const World& world, Vector3 start1, Vector3 end1, const std::vector<ShellID>& ignored_shells) :
            world(world), start1(start1), end1(end1), ignored_shells(ignored_shells), is_limited(false), max_radius(1e6) {}

    void operator()(typename TCollector::const_iterator i, const Vector3&)
    {
        auto shellID = (*i).first;
        auto shell = (*i).second;

        if(std::find(ignored_shells.begin(),ignored_shells.end(), shellID) != ignored_shells.end())
        {
            // Skip shells in ignored list
            return;
        }

        Vector3 start2, end2;
        double distance = 0.0, other_radius = 0.0;
        if(shell.shape() == Shell::Shape::SPHERE)
        {
            auto sphere = shell.get_sphere();
            other_radius = sphere.radius();

            // Place the cylinder closest to the sphere, given a cyclic world.
            start1 = world.cyclic_transpose(start1, sphere.position());
            end1 = world.cyclic_transpose(end1, sphere.position());

            distance = std::sqrt(squared_distance::line_segment_point(start1, end1, sphere.position()));

            if (feq(distance, 0.0, 1e-12))
            {
                return;
            }
        }
        else if (shell.shape() == Shell::Shape::CYLINDER)
        {
            auto cylinder = shell.get_cylinder();
            start2 = cylinder.position() - cylinder.half_length() * cylinder.unit_z();
            end2 = cylinder.position() + cylinder.half_length() * cylinder.unit_z();
            other_radius = cylinder.radius();

            // Place the cylinder closest to the other cylinder, given a cyclic world.
            start1 = world.cyclic_transpose(start1, cylinder.position());
            end1 = world.cyclic_transpose(end1, cylinder.position());

            distance = std::sqrt(squared_distance::line_segments(start1, end1, start2, end2));

            if (feq(distance, 0.0, 1e-12))
            {
                return;
            }
        }
        else
        {
            THROW_EXCEPTION(illegal_argument, "Shell shape to scale to is neither spherical or cylindrical")
        }

        auto available = distance - other_radius;
        if(available < max_radius)
        {
            max_radius = available;
            is_limited = true;
            limiting_shell_type = shell.code();
        }
    }
    const World& world;
    Vector3 start1;
    Vector3 end1;
    const std::vector<ShellID>& ignored_shells;
    double max_radius;
    bool is_limited; // true if the radius is limited by a neighbouring shell
    Shell::Code limiting_shell_type; // shell type of neighbouring cell that limits radius
};

/**
 * Finds a maximal dynamical cylinder height (where total cylinder height is static height + dynamical height)
 * for a given cylinder, such that it does not touch any neighbouring structures or domains.
 *
 * In case of a Mixed2D3DPair or SinglePlanarInteraction, the static height is the distance from the 3D-diffusing
 * particle to their PlanarSurface + the length by which the cylinder exceeds beyond the other side of the PlanarSurface,
 * and the dynamic height is the length from the 3D-diffusing particle to the cylinder cap away from the PlanarSurface.
 *
 * @tparam TMatrixSpace
 * @param base_pos Position of the bottom of the cylinder. That is, the side of the 2D-diffusing particle.
 * @param unit_z Unit vector pointing from the bottom to the top of the cylinder.
 * @param height_static The static length of the cylinder, as described above.
 * @param scaling_factor The factor such that height_dynamic = radius / scaling_factor, that ensures isotropic diffusion of the IV.
 * @param matrixSpace The matrix space that divides the total world region into cells.
 * @param structures All structures in the world region.
 * @param ignored_structures A list of structure IDs that will not be compared against (e.g. attached PlanarSurfaces, and the world).
 * @return The maximal dynamical height of the cylinder, such that no neighbouring structures or domains are touched.
 */
template <typename TMatrixSpace>
double scaling::find_maximal_cylinder_height(Vector3 base_pos, Vector3 unit_z, double height_static,
                                             double scaling_factor, const TMatrixSpace& matrixSpace,
                                             const World& world,
                                             const std::vector<StructureID>& ignored_structures,
                                             const std::vector<ShellID>& ignored_shells,
                                             bool social_scaling)
{
    double max_dynamic_height = 1e6;
    auto structures = world.get_structures();

    for(const auto& structure : structures)
    {
        if (std::find(ignored_structures.begin(), ignored_structures.end(), structure->id()) == ignored_structures.end()) {
            auto plane = dynamic_cast<PlanarSurface *>(structure.get());
            if (plane == nullptr) {
                THROW_EXCEPTION(not_implemented,
                                "Scaling cylinders up to structures other than PlanarSurfaces is not yet implemented.");
            }

            auto transposed_base = world.cyclic_transpose(base_pos, plane->position());
            auto calculated_max_height = find_maximal_cylinder_height_to_plane(base_pos, unit_z,
                                                                               height_static,
                                                                               scaling_factor, plane);
            max_dynamic_height = std::min(max_dynamic_height, calculated_max_height);
        }
    }

    max_height_collector<TMatrixSpace> col(world, base_pos, unit_z, height_static, scaling_factor, ignored_shells);
    matrixSpace.each_neighbor_cyclic(matrixSpace.index(base_pos), col);
    max_dynamic_height = std::min(max_dynamic_height, col.max_height);


    if(social_scaling && max_dynamic_height != 1e6)
    {
        return social_correction(max_dynamic_height, col.limiting_shell_type);
    }

    return max_dynamic_height;
}


template <typename TMatrixSpace>
double scaling::find_maximal_cylinder_radius(Vector3 start1, Vector3 end1, const TMatrixSpace& matrixSpace,
                                             const World& world,
                                             const std::vector<StructureID>& ignored_structures,
                                             const std::vector<ShellID>& ignored_shells,
                                             bool social_scaling) {
    double max_radius = 1e6;
    auto structures = world.get_structures();

    for (const auto &structure : structures) {
        if (std::find(ignored_structures.begin(), ignored_structures.end(), structure->id()) ==
            ignored_structures.end()) {
            auto plane = dynamic_cast<PlanarSurface *>(structure.get());
            if (plane == nullptr) {
                THROW_EXCEPTION(not_implemented,
                                "Scaling cylinders up to structures other than PlanarSurfaces is not yet implemented.");
            }

            auto transposed_start1 = world.cyclic_transpose(start1, plane->position());
            auto transposed_end1 = world.cyclic_transpose(end1, plane->position());

            auto start2 = plane->project_point(transposed_start1).first;
            auto end2 = plane->project_point(transposed_end1).first;
            auto distance = std::sqrt(squared_distance::line_segments(transposed_start1, transposed_end1, start2, end2));
            max_radius = std::min(max_radius, distance);
        }
    }

    max_radius_collector<TMatrixSpace> col(world, start1, end1, ignored_shells);
    matrixSpace.each_neighbor_cyclic(matrixSpace.index(start1), col);
    max_radius = std::min(max_radius, col.max_radius);

    if(social_scaling && max_radius < 1e6)
    {
        return social_correction(max_radius, col.limiting_shell_type);
    }

    return max_radius;
}