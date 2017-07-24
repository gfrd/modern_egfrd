#ifndef REGION_HPP
#define REGION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "Structure.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

/* NOTE: the class Region is exactly identical to the Surface class (so why need this?)


template<typename Tshape_>
class Region : public Structure           // ParticleSimulationStructure<Ttraits_> <- removed class in hierarchy but may be needed to stuff it back in again
{
public:
typedef Tshape_                                   shape_type;
typedef std::pair<Vector3, bool>                  position_flag_pair_type;
typedef std::pair<double, double>                 DoublePair;
typedef std::pair<Vector3, DoublePair>  projected_type;

Region(std::string const& name, SpeciesTypeID const& sid, StructureID const& parent_struct_id, shape_type const& shape)
: Structure(name, sid, parent_struct_id), shape_(shape) {}

virtual ~Region() {}

shape_type const& shape() const { return shape_; }

virtual bool operator==(Structure const& rhs) const override
{
Region const* _rhs(dynamic_cast<Region const*>(&rhs));
return _rhs && Structure::operator==(rhs) && shape_ == _rhs->shape_;
}

virtual projected_type project_point(Vector3 const& pos) const override
{
return shape_.project_point(pos);
}

virtual projected_type project_point_on_surface(Vector3 const& pos) const override
{
return shape_.project_point_on_surface(pos);
}

virtual double distance(Vector3 const& pos) const override
{
return shape_.distance(pos);
}

virtual Vector3 const& position() const override
{
return shape_.position();
}

virtual position_flag_pair_type deflect(Vector3 const& pos0, Vector3 const& displacement) const override
{
return shape_.deflect(pos0, displacement);
}

//virtual Vector3 deflect_back(Vector3 const& pos, Vector3 const& u_z) const
//{
//   return shape_.deflect_back(pos, u_z);
//}

virtual std::size_t hash() const override
{
return Structure::hash() ^ std::hash<shape_type>()();
}

virtual std::string as_string() const override
{
std::ostringstream out;
out << "Region(" << domainID_ << ", " << sid_ << ", " << name_ << " : " + shape_ + ")";
return out.str();
}

protected:
shape_type shape_;
};
*/

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* REGION_HPP */
