#include "PlanarSurface.hpp"
#include "StructureContainer.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

PlanarSurface::position_structid_pair PlanarSurface::apply_boundary(const position_structid_pair& pos_structure_id, const StructureContainer& sc) const
{
    Vector2 half_extents = shape().half_extent();
    Vector3 pos_vector = pos_structure_id.first - shape().position();

    double   cx = Vector3::dot(pos_vector, shape().unit_x());
    double   cy = Vector3::dot(pos_vector, shape().unit_y());
    //double   cz = Vector3::dot(pos_vector, shape().unit_z());

    //// Check for (currently unsupported) self-connections
    //for (int i = 0; i < 4; i++)
    //{
    //    auto neighbor_id_vector = sc.get_neighbor_info(*this, 0);
    //    if (neighbor_id_vector.first == pos_structure_id.second)
    //    {
    //        //log_.warn("Plane is connected to itself, which is unsupported; no boundary condition will be applied at this call.");
    //        return pos_structure_id;
    //    }
    //};

    if ((std::abs(cx) <= half_extents.X() && std::abs(cy) <= half_extents.Y()))
        return pos_structure_id; // we are still in the plane (did not pass any of the boundaries)

    StructureID new_id = pos_structure_id.second;   // initialized with old structure ID
    Vector3 new_pos = position();

    if (half_extents.Y() * std::abs(cx) < half_extents.X() * std::abs(cy))
    {
        // Top or bottom (and may also be in one of the corners)
        if (half_extents.Y() < cy)
        {
            // we are at the 'top' of the plane (side nr. 0)
            // the 'vector' is the unit vector pointing from the edge between the two planes to the center of the plane
            auto neighbor_id_vector = sc.get_neighbor_info(*this, 0);
            new_id = neighbor_id_vector.first;

            if (neighbor_id_vector.second != Vector3::null)
            {
                Vector3 neighbor_plane_par = shape().unit_x() * cx;
                Vector3 neighbor_plane_inl = shape().unit_y() * half_extents.Y() + neighbor_id_vector.second * (cy - half_extents.Y());
                new_pos = position() + neighbor_plane_par + neighbor_plane_inl;
            }
        }
        else if (cy < -half_extents.Y())
        {
            // we are at the 'bottom' of the plane (side nr. 1)
            auto neighbor_id_vector = sc.get_neighbor_info(*this, 1);
            new_id = neighbor_id_vector.first;

            if (neighbor_id_vector.second != Vector3::null)
            {
                Vector3 neighbor_plane_par = shape().unit_x() * cx;
                Vector3 neighbor_plane_inl = shape().unit_y() * -half_extents.Y() + neighbor_id_vector.second * (-cy - half_extents.Y());
                new_pos = position() + neighbor_plane_par + neighbor_plane_inl;
            }
        }
    }
    else // half_extents.Y()*abs(cx) >= half_extents[0]*abs(cy)
    {
        if (cx < -half_extents.X())
        {
            // we are at the 'left' of the plane (side nr. 2)
            auto neighbor_id_vector = sc.get_neighbor_info(*this, 2);
            new_id = neighbor_id_vector.first;

            if (neighbor_id_vector.second != Vector3::null)
            {
                Vector3 neighbor_plane_par = shape().unit_y() * cy;
                Vector3 neighbor_plane_inl = shape().unit_x() * -half_extents.X() + neighbor_id_vector.second * (-cx - half_extents.X());
                new_pos = position() + neighbor_plane_par + neighbor_plane_inl;
            }
        }
        else if (half_extents.X() < cx)
        {
            // we are at the 'right' of the plane (side nr. 3)
            auto neighbor_id_vector = sc.get_neighbor_info(*this, 3);
            new_id = neighbor_id_vector.first;

            if (neighbor_id_vector.second != Vector3::null)
            {
                Vector3 neighbor_plane_par = shape().unit_y() *  cy;
                Vector3 neighbor_plane_inl = shape().unit_x() * half_extents.X() + neighbor_id_vector.second * (cx - half_extents.X());
                new_pos = position() + neighbor_plane_par + neighbor_plane_inl;
            }
        }
    }

    // Check if we are in one of the corners. If yes -> do another round of border crossing from neighboring plane.
    // Only do this if the particle does not end up on another side of the same plane (the scenario of one plane 
    // with periodic BCs), as this may cause infinite loops.
    if (new_id != pos_structure_id.second && std::abs(cx) > half_extents.X() && std::abs(cy) > half_extents.Y())
        return sc.get_structure(new_id)->apply_boundary(std::make_pair(new_pos, new_id), sc);

    return std::make_pair(new_pos, new_id);
}

// --------------------------------------------------------------------------------------------------------------------------------

PlanarSurface::position_structid_pair PlanarSurface::cyclic_transpose(const position_structid_pair& pos_structure_id, const StructureContainer& sc) const
{
    const PlanarSurface& planar_surface = *this;
    typedef PlanarSurface::shape_type    plane_type;
    typedef std::pair<StructureID, Vector3>     neighbor_id_vector;

    if (pos_structure_id.second == planar_surface.id())
        return pos_structure_id;

    Vector3 new_pos_par;
    Vector3 new_pos_inl;

    int side_nr;
    neighbor_id_vector niv;
    for (side_nr = 0; side_nr < 4; ++side_nr)
    {
        niv = sc.get_neighbor_info(planar_surface, side_nr);
        if (niv.first == pos_structure_id.second)
            // we have found the correct side
            break;
    }
    // TODO use iterator so that we can check whether match was found or that we escaped by reaching the end of the list.

    const plane_type plane(planar_surface.shape());
    Vector2 half_extents(plane.half_extent());
    switch (side_nr)
    {
    case 0:     // top side
    {
        const Vector3 pos((pos_structure_id.first - (plane.position() + plane.unit_y() * half_extents.Y())));
        const double plane_par(Vector3::dot(pos, plane.unit_x()));     // TODO replace by subtraction?
        const double plane_inl(Vector3::dot(pos, niv.second));
        new_pos_par = plane.unit_x()* plane_par;
        new_pos_inl = plane.unit_y()* (half_extents.Y() + plane_inl);
        break;
    }
    case 1:     // bottom
    {
        const Vector3 pos((pos_structure_id.first - (plane.position() + plane.unit_y()* -half_extents.Y())));
        const double plane_par(Vector3::dot(pos, plane.unit_x()));     // TODO replace by subtraction?
        const double plane_inl(Vector3::dot(pos, niv.second));
        new_pos_par = plane.unit_x()* plane_par;
        new_pos_inl = plane.unit_y()* -(half_extents.Y() + plane_inl);
        break;
    }
    case 2:     // left
    {
        const Vector3 pos((pos_structure_id.first - (plane.position() + plane.unit_x()* -half_extents.X())));
        const double plane_par(Vector3::dot(pos, plane.unit_y()));     // TODO replace by subtraction?
        const double plane_inl(Vector3::dot(pos, niv.second));
        new_pos_par = plane.unit_y()* plane_par;
        new_pos_inl = plane.unit_x()* -(half_extents.X() + plane_inl);
        break;
    }
    case 3:     // right
    {
        const Vector3 pos((pos_structure_id.first - (plane.position() + plane.unit_x()* half_extents.X())));
        const double plane_par(Vector3::dot(pos, plane.unit_y()));     // TODO replace by subtraction?
        const double plane_inl(Vector3::dot(pos, niv.second));
        new_pos_par = plane.unit_y()* plane_par;
        new_pos_inl = plane.unit_x()* (half_extents.X() + plane_inl);
        break;
    }
    }
    const Vector3 new_pos((plane.position() + (new_pos_par + new_pos_inl)));
    return std::make_pair(new_pos, planar_surface.id());
}

// --------------------------------------------------------------------------------------------------------------------------------

