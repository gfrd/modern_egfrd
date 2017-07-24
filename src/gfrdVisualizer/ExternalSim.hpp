#ifndef EXTERNALSIM_HPP
#define EXTERNALSIM_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <vector>
#include <cstring>
#include <Particle.hpp>
#include <Single.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

class ExtSim
{
public:
   ExtSim() = default;

   bool active() const { return std::strlen(szFilename) != 0; }
   double time() const { return time_; }
   size_t num_steps() const { return num_steps_; }
   Vector3 world_size() const { return Vector3(world_size_, world_size_, world_size_); }
   double cell_size() const { return world_size_ / matrix_size_; }
   std::array<uint, 3> matrix_size() const { return std::array<uint, 3>({ matrix_size_,matrix_size_,matrix_size_ }); }
   const std::vector<Particle::particle_id_pair>& get_particles() const { return particles_; }

   int Section() const noexcept { return section + 1; }
   int Sections()const noexcept { return sections; }

   StructureContainer::structures_range get_structures() const { return gi::iteratorRange<StructureContainer::structure_map, std::shared_ptr<Structure>, gi::SelectSecond<std::shared_ptr<Structure>>>(structures_); }
   const std::vector<std::pair<Shell, int>>& get_domains() const { return domains_; }
   const char* get_filename() const { return szFilename; }

   bool readSimFile(const char*filename);
   bool refresh();
   int SelectSection(int selSection);

private:
   std::vector<Particle::particle_id_pair> particles_;
   std::vector<std::pair<Shell,int>> domains_;
   std::map<StructureID, std::shared_ptr<Structure>> structures_;
   double time_;
   size_t num_steps_;
   double world_size_;
   uint matrix_size_;
   char szFilename[1024];
#if defined(_MSC_VER)
   uint last_write_hi_;
   uint last_write_lo_;
#else
#endif
   uint sections;
   int section;
};

// --------------------------------------------------------------------------------------------------------------------------------

extern ExtSim extSim;

// --------------------------------------------------------------------------------------------------------------------------------

#endif