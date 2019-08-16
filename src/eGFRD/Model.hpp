#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <algorithm>
#include <map>
#include "genericIterator.hpp"
#include "SerialIDGenerator.hpp"
#include "SpeciesTypeID.hpp"
#include "StructureTypeID.hpp"
#include "SpeciesType.hpp"
#include "makeString.hpp"
#include "StructureType.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Model
{
   using species_map = std::map<SpeciesTypeID, SpeciesType>;
   using structure_map = std::map<StructureTypeID, StructureType>;
   using SpeciesIteratorRange = gi::iteratorRange<species_map>;
   using StructureIteratorRange = gi::iteratorRange<structure_map>;

public:
   Model() noexcept
   {
      auto default_structure_type = StructureType("world");
      add_structure_type(default_structure_type);
      default_structure_type_id_ = default_structure_type.id();
   }

   template <typename speciesType>     // universal reference, accept l-values and r-values
   SpeciesTypeID add_species_type(speciesType&& species)
   {
      auto id = spgen_();
      species.set_id(id);
      species_map_.insert(std::make_pair(id, std::forward<SpeciesType>(species)));
      return id;
   }

   const SpeciesType& get_species_type_by_id(const SpeciesTypeID id) const
   {
      auto i = species_map_.find(id);
      if (species_map_.end() == i) throw not_found(make_string() << "SpeciesType with id=" << id << " not found!");
      return i->second;
   }

   SpeciesTypeID get_species_type_id_by_name(const std::string name) const
   {
      auto i = std::find_if(species_map_.begin(), species_map_.end(), [name](const species_map::value_type& n) { return n.second.name() == name; });
      if (species_map_.end() == i) throw not_found(make_string() << "SpeciesType with name='" << name << "' not found!");
      return i->first;
   }

   SpeciesIteratorRange get_species() const { return SpeciesIteratorRange(species_map_); }

   StructureTypeID add_structure_type(StructureType&& structure)
   {
      structure.set_id(stgen_());
      structure_map_.insert(std::make_pair(structure.id(), std::move(structure)));
      return structure.id();
   }

   void add_structure_type(StructureType& structure)
   {
      structure.set_id(stgen_());
      structure_map_.insert(std::make_pair(structure.id(), structure));
   }

   const StructureType& get_structure_type_by_id(const StructureTypeID id) const
   {
      structure_map::const_iterator i(structure_map_.find(id));
      if (structure_map_.end() == i) throw not_found(make_string() << "StructureType with id=" << id << "not found!");
      return (*i).second;
   }

    StructureTypeID get_structure_type_id_by_name(const std::string name) const
    {
        auto i = std::find_if(structure_map_.begin(), structure_map_.end(), [name](const structure_map::value_type& n) { return n.second.name() == name; });
        if (structure_map_.end() == i) { return StructureTypeID(0); }
        return i->first;
    }

   StructureTypeID get_def_structure_type_id() const { return default_structure_type_id_; }

   StructureIteratorRange structures() const { return StructureIteratorRange(structure_map_); }

private:
   SerialIDGenerator<SpeciesTypeID>    spgen_;                             // The id generator which makes sure that all species-/structureTypes have a unique id
   SerialIDGenerator<StructureTypeID>  stgen_;                             // The id generator which makes sure that all species-/structureTypes have a unique id
   species_map                         species_map_;                       // mapping: SpeciesTypeID->SpeciesType
   structure_map                       structure_map_;                     // mapping: StructureTypeID->StructureType
   StructureTypeID                     default_structure_type_id_;         // The id of the default structure_type ("world")
};

// --------------------------------------------------------------------------------------------------------------------------------
