#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "Vector3.hpp"
#include "StructureID.hpp"
#include "randomNumberGenerator.hpp"
#include "StructureTypeID.hpp"
#include "makeString.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct SpeciesTypeID;
class SpeciesType;
// Forward declarations
class StructureContainer;
class CuboidalRegion;
class CylindricalSurface;
class SphericalSurface;
class DiskSurface;
class PlanarSurface;

// --------------------------------------------------------------------------------------------------------------------------------

class Structure
{
public:
   typedef std::pair<double, double>                                             DoublePair;
   typedef std::pair<Vector3, DoublePair>                                        projected_type;
   typedef std::pair<Vector3, Vector3>                                           position_pair_type;
   typedef std::pair<Vector3, bool>                                              position_flag_pair_type;
   typedef std::pair<Vector3, StructureID>                                       position_structid_pair;
   typedef std::pair<position_structid_pair, position_structid_pair>             position_structid_pair_pair;

   Structure(const std::string& name, const StructureTypeID sid, const StructureID pid) noexcept
      : name_(name), sid_(sid), parent_id_(pid) {}

   virtual ~Structure() = default;

   StructureID id() const { return id_; }

   const std::string& name() const { return name_; }

   StructureTypeID sid() const { return sid_; }

   StructureID parent_id() const { return parent_id_; }

   virtual bool operator==(const Structure& rhs) const
   {
      return id_ == rhs.id() && sid_ == rhs.sid();
   }

   bool operator!=(const Structure& rhs) const
   {
      return !operator==(rhs);
   }

   virtual Vector3 random_position(RandomNumberGenerator& rng) const = 0;
   virtual Vector3 random_vector(double r, RandomNumberGenerator& rng) const = 0;

   // Methods used in the 'new' BDPropagator
   virtual Vector3 bd_displacement(double mean, double r, RandomNumberGenerator& rng) const = 0;
   virtual double newBD_distance(const Vector3& new_pos, double radius, const Vector3& old_pos, double sigma) const = 0;

   // TODO this are just functions->move somewhere else
//   virtual double get_1D_rate_geminate(double k, double r01) const = 0;
   virtual double get_1D_rate_surface(double k, double r0) const = 0;
   virtual double particle_reaction_volume(double r01, double rl) const = 0;
   virtual double surface_reaction_volume(double r0, double rl) const = 0;     // This does contain a surface dependent component.
//
//   // Methods used to calculate dissociation positions
//   virtual Vector3 surface_dissociation_vector(RandomNumberGenerator& rng, double r0, double rl) const = 0;
//   virtual Vector3 surface_dissociation_unit_vector(RandomNumberGenerator& rng) const = 0;
////   virtual position_pair_type geminate_dissociation_positions(RandomNumberGenerator& rng, const SpeciesType& s0, const SpeciesType& s1, const Vector3& op, double rl) const = 0;
////   virtual position_pair_type special_geminate_dissociation_positions(RandomNumberGenerator& rng, const SpeciesType& s_surf, const SpeciesType& s_bulk, const Vector3& op_surf, double rl) const = 0;
//
//   // General method for getting some measures/info
   virtual projected_type project_point(const Vector3& pos) const = 0;
   //   virtual projected_type project_point_on_surface(const Vector3& pos) const = 0;
   virtual double distance(const Vector3& pos) const = 0;
   virtual const Vector3& position() const = 0;
   virtual Vector3 side_comparison_vector() const = 0;
   //
   //   // Methods used for edge crossing (only for the planes so far)
   //   virtual position_flag_pair_type deflect(const Vector3& pos0, const Vector3& displacement) const = 0;
      //    virtual Vector3 deflect_back(const Vector3& pos, const Vector3& u_z) const = 0;

   virtual position_structid_pair apply_boundary(const position_structid_pair& pos_struct_id, const StructureContainer& structure_container) const = 0;
   virtual position_structid_pair cyclic_transpose(const position_structid_pair& pos_struct_id, const StructureContainer& structure_container) const = 0;

   // *** Structure functions dynamic dispatch ***
   // 
   // FIXME For now the second dispatch requires the helper functions to be defined here for each structure type separately.
   // This is because C++ does not allow virtual templates. The current solution is functional, but ugly, and may be replaced
   // by a more elegant solution in the future.
   // 
   // *** 1 *** - Producing one new position
   // First dispatch
   virtual position_structid_pair get_pos_sid_pair(const Structure& target_structure, const Vector3& position, double offset, double rl, RandomNumberGenerator& rng) const = 0;
   // Second dispatch
   // The helper functions are the dispatch acceptors and have to be declared for each derived structure class because C++ does not support virtual templates.
   //virtual position_structid_pair get_pos_sid_pair_helper(CuboidalRegion const& origin_structure, const Vector3& position,double offset, double rl, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair get_pos_sid_pair_helper(SphericalSurface const& origin_structure, const Vector3& position,double offset, double rl, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair get_pos_sid_pair_helper(CylindricalSurface const& origin_structure, const Vector3& position,double offset, double rl, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair get_pos_sid_pair_helper(DiskSurface const& origin_structure, const Vector3& position, double offset, double rl, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair get_pos_sid_pair_helper(PlanarSurface const& origin_structure, const Vector3& position,double offset, double rl, RandomNumberGenerator& rng) const = 0;
   //// *** 2 *** - Producing two new positions
   //// First dispatch
   virtual position_structid_pair_pair get_pos_sid_pair_pair(const Structure& target_structure, const Vector3& position, const SpeciesType& s_orig, const SpeciesType& s_targ, double rl, RandomNumberGenerator& rng) const = 0;
   //// Second dispatch
   //// The helper functions are dispatch acceptors and have to be declared for each derived structure class because C++ does not support virtual templates.
   //virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(CuboidalRegion const& origin_structure, const Vector3& position, const SpeciesType& s_orig, const SpeciesType& s_targ, double rl, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(SphericalSurface const& origin_structure, const Vector3& position, const SpeciesType& s_orig, const SpeciesType& s_targ, double rl, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(CylindricalSurface const& origin_structure, const Vector3& position, const SpeciesType& s_orig, const SpeciesType& s_targ, double rl, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(DiskSurface const& origin_structure, const Vector3& position, const SpeciesType& s_orig, const SpeciesType& s_targ, double rl, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair_pair get_pos_sid_pair_pair_helper(PlanarSurface const& origin_structure, const Vector3& position, const SpeciesType& s_orig, const SpeciesType& s_targ, double rl, RandomNumberGenerator& rng) const = 0;

   // *** 3 *** - Pair reactions => two origin structures
   // The following functions handle the case of two origin structures.
   // First again the (C++) structure types have to be determined by a double dispatch.
   // As a next step, the helper function has to determine which of the two structures
   // is the target structure.
   // For now, the target structure is either:
   //   - the lower hierarchy level structure, i.e. one of the origin structures has
   //     to be the daughter structure of the other and the particle will end up on
   //     the daughter structure; or:
   //   - in case of equal structure type id's it can end up on either origin structure
   //     and apply_boundary will handle the right placement afterwards.

   // First dispatch
   // This is called as a method of origin_structure1 with origin_structure2 as an argument.
   virtual position_structid_pair get_pos_sid_pair_2o(const Structure& origin_structure2, const SpeciesTypeID target_sid, const Vector3& CoM, double offset, double reaction_length, RandomNumberGenerator& rng) const = 0;
   //     // Some convenient method overloading; this is just a redirect to the above                                       
   //     virtual position_structid_pair get_pos_sid_pair(const Structure& origin_structure2, const SpeciesTypeID& target_sid, const Vector3& CoM,
   //                                                          double offset, double reaction_length, RandomNumberGenerator& rng) const = 0;
   // Second dispatch
   // The helper functions are the dispatch acceptors and have to be declared for each derived structure class because C++ does not support virtual templates.
   //virtual position_structid_pair get_pos_sid_pair_2o_helper(CuboidalRegion const& origin_structure1, const StructureID& target_sid, const Vector3& CoM, double offset, double reaction_length, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair get_pos_sid_pair_2o_helper(SphericalSurface const& origin_structure1, const SpeciesTypeID& target_sid, const Vector3& CoM,double offset, double reaction_length, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair get_pos_sid_pair_2o_helper(CylindricalSurface const& origin_structure1, const SpeciesTypeID& target_sid, const Vector3& CoM, double offset, double reaction_length, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair get_pos_sid_pair_2o_helper(DiskSurface const& origin_structure1, const SpeciesTypeID& target_sid, const Vector3& CoM, double offset, double reaction_length, RandomNumberGenerator& rng) const = 0;
   //virtual position_structid_pair get_pos_sid_pair_2o_helper(PlanarSurface const& origin_structure1, const SpeciesTypeID& target_sid, const Vector3& CoM, double offset, double reaction_length, RandomNumberGenerator& rng) const = 0;

   // Some further helper functions used by template<typename Tstruct_> get_pos_sid_pair_helper_two_origins_any(...),
   // which is the final dispatch template defined in each of the derived classes and makes use of the two following checker functions:
   bool is_parent_of_or_has_same_sid_as(const Structure& s) const
   {
      return s.parent_id() == id_ || s.sid() == sid_;
   }

   bool has_valid_target_sid(const StructureTypeID target_sid) const
   {
      return sid_ == target_sid;
   }

   //     // TODO
   //     // *** 4 *** - Generalized functions for pair reactions => two origin structures and one target_structure
   //     // This introduces a triple dynamic dispatch, overloading method call structure.get_pos_sid_pair once more.
   //     // NOTE: As yet these methods are unused but might prove useful in the future.
   //     virtual position_structid_pair get_pos_sid_pair(const Structure& origin_structure2, const Structure& target_structure, const Vector3& position,
   //                                                          double offset, double reaction_length, RandomNumberGenerator& rng) const = 0;
   //     template <typename Tstruct1_>
   //     position_structid_pair get_pos_sid_pair_helper1(Tstruct1_ const& origin_structure1, const Structure& target_structure, const Vector3& position,
   //                                                          double offset, double reaction_length, RandomNumberGenerator& rng) const
   //     {
   //         return target_structure.get_pos_sid_pair_helper2(origin_structure1, *this, position, offset, reaction_length, rng);
   //     }
   //     template <typename Tstruct1_, typename Tstruct2_>
   //     position_structid_pair get_pos_sid_pair_helper2(Tstruct1_ const& origin_structure1, Tstruct2_ const& origin_structure2, const Vector3& position,
   //                                                          double offset, double reaction_length, RandomNumberGenerator& rng) const
   //     {
   //         StructureID    this_id( this->id );
   //         StructureID    os1_id( origin_structure1.id );
   //         StructureID    os2_id( origin_structure2.id );
   //
   //         if(os1_id == this_id)
   //             // Dispatch to function with well-defined typing
   //             return ::get_pos_sid_pair(origin_structure2, *this, position, offset, reaction_length, rng);
   //         
   //         else if(os2_id == this_id)
   //             // Dispatch to function with well-defined typing
   //             return ::get_pos_sid_pair(origin_structure1, *this, position, offset, reaction_length, rng);
   //         
   //         else
   //             throw propagation_error("Target structure must be one of the origin structures for pair reaction.");
   //     }
   //     
   //     NOTE: The template based variant will not work! The helper methods have to be defined for each structure type separately or in a smarter way!


   virtual std::size_t hash() const
   {
      return std::hash<std::string>()(name_) ^ std::hash<StructureTypeID>()(sid_);
   }

   virtual std::string as_string() const = 0;      // runtime parameters (name, id's position, etc)

   virtual const char* type_name() const { return "Structure"; };     // type implementation name (Structure, PlaneSurface, CuboidalRegion, etc)


protected:
   friend class StructureContainer;    // only he can set the id!
   void set_id(const StructureID id) { id_ = id; }

   friend class Persistence;

   std::string       name_;
   StructureTypeID   sid_;
   StructureID       id_;
   StructureID       parent_id_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const Structure& s)
{
   stream << s.type_name() << "{" "'" << s.name() << "', " << s.id() << ", parent=" << s.parent_id() << ", " << s.sid() << s.as_string() << "}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

namespace std {
   template<>
   struct hash < Structure >
   {
      std::size_t operator()(const Structure& s) { return s.hash(); }
   };
} // namespace std

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* STRUCTURE_HPP */
