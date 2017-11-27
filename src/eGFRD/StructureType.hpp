#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <ostream>
#include "exceptions.hpp"
#include "StructureTypeID.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class StructureType
{
public:
   StructureType() = default;                                   // use default construct
   StructureType(const StructureType&) = default;               // copy operation
   StructureType& operator=(const StructureType&) = default;
   StructureType(StructureType&&) = default;                    // move operation
   StructureType& operator=(StructureType&&) = default;

   StructureType(const std::string& name) noexcept : name_(name) {}

   StructureTypeID id() const
   {
      if (!id_) throw illegal_state("Not bound to Model");
      return id_;
   };

   // Get the name
   const std::string & name() const { return name_; }

protected:
   friend class Model;  // class that can set id
   void set_id(const StructureTypeID id) { id_ = id; }

   friend class Persistence;

private:
   StructureTypeID id_;
   std::string name_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const StructureType& s)
{
   stream << "StructureType{" << s.id() << " '" << s.name() << "'}";
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------
