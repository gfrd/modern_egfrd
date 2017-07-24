
#include <fstream>
#include <algorithm>
#include "ExternalSim.hpp"

#if defined(_MSC_VER)
#define NOMINMAX
#include "windows.h"
#endif

// --------------------------------------------------------------------------------------------------------------------------------

ExtSim extSim;

// --------------------------------------------------------------------------------------------------------------------------------

bool ExtSim::readSimFile(const char*filename)
{
   std::ifstream infile(filename);
   std::string line;

   std::vector<Particle::particle_id_pair> particles;
   std::vector<std::pair<Shell, int>> domains;
   std::map<StructureID, std::shared_ptr<Structure>> structures;
   double time = 0.0;
   size_t num_steps = 0;
   double world_size = 0.0;
   uint matrix_size = 0;

   try
   {
      std::ifstream::pos_type last_pos(0), start_section;
      int section_count = 0;
      while (std::getline(infile, line))
      {
         size_t pos = line.find("sim_step");
         if (pos == 0 && section == section_count) start_section = last_pos;    // capture pos of selected section (or last section when -1)
         last_pos = infile.tellg();

         pos = line.find("End");
         if (pos == 0) section_count++;
      }

      infile.close();
      if (section_count == 0) { sections = 0; return false; }

      infile.open(filename);
      infile.seekg(start_section, std::ios_base::beg);
      last_pos = infile.tellg();

      while (std::getline(infile, line))
      {
         size_t pos = line.find("sim_step");
         if (pos == 0)
            sscanf(line.c_str(), "sim_step=%zu\tsim_time=%le\tworld_size=%le\tmatrix_size=%d", &num_steps, &time, &world_size, &matrix_size);

         pos = line.find("Particle");
         if (pos == 0)
         {
            int id, sid; double r, x, y, z;
            sscanf(line.c_str(), "Particle=%d\tRadius=%le\tPosition=%le,%le,%le\tSid=%d", &id, &r, &x, &y, &z, &sid);
            particles.emplace_back(std::make_pair(ParticleID(id), Particle(SpeciesTypeID(sid), Sphere(Vector3(x, y, z), r), StructureID(0), 0.0, 0.0)));
         }

         pos = line.find("Structure");
         if (pos == 0)
         {
            int id, sid; char name[64], shape[16]; double x, y, z;
            sscanf(line.c_str(), "Structure=%d\tName=%s\tSid=%d\tShape=%s\tPosition=%le,%le,%le", &id, name, &sid, shape, &x, &y, &z);
            if (strcmp(shape, "CUBE") == 0)
            {
               double hx, hy, hz;
               sscanf(line.c_str(), "%*s\t%*s\t%*s\t%*s\t%*s\tHalfExtent=%le,%le,%le", &hx, &hy, &hz);
               structures[StructureID(id)] = (std::make_unique<CuboidalRegion>(std::string(name), StructureTypeID(sid), StructureID(id), Box(Vector3(x, y, z), Vector3(hx, hy, hz))));
            }
            else if (strcmp(shape, "PLANE") == 0)
            {
               double hx, hy, uxx, uxy, uxz, uyx, uyy, uyz; int o;
               sscanf(line.c_str(), "%*s\t%*s\t%*s\t%*s\t%*s\tHalfExtent=%le,%le\tUnitX=%le,%le,%le\tUnitY=%le,%le,%le\tOneSide=%d", &hx, &hy, &uxx, &uxy, &uxz, &uyx, &uyy, &uyz, &o);
               structures[StructureID(id)] = (std::make_unique<PlanarSurface>(std::string(name), StructureTypeID(sid), StructureID(id), Plane(Vector3(x, y, z), Vector3(uxx, uxy, uxz), Vector3(uyx, uyy, uyz), hx, hy, o != 0)));
            }
            else if (strcmp(shape, "SPHERE") == 0)
            {
               double r;
               sscanf(line.c_str(), "%*s\t%*s\t%*s\t%*s\t%*s\tRadius=%le", &r);
               structures[StructureID(id)] = (std::make_unique<SphericalSurface>(std::string(name), StructureTypeID(sid), StructureID(id), Sphere(Vector3(x, y, z), r)));
            }
            else if (strcmp(shape, "CYLINDER") == 0)
            {
               double r, ux, uy, uz, hl;
               sscanf(line.c_str(), "%*s\t%*s\t%*s\t%*s\t%*s\tRadius=%le\tUnitZ=%le,%le,%le\tHalfLength=%le", &r, &ux, &uy, &uz, &hl);
               structures[StructureID(id)] = (std::make_unique<CylindricalSurface>(std::string(name), StructureTypeID(sid), StructureID(id), Cylinder(Vector3(x, y, z), r, Vector3(ux, uy, uz), hl)));
            }
            else if (strcmp(shape, "DISK") == 0)
            {
               double r, ux, uy, uz;
               sscanf(line.c_str(), "%*s\t%*s\t%*s\t%*s\t%*s\tRadius=%le\tUnitZ=%le,%le,%le", &r, &ux, &uy, &uz);
               structures[StructureID(id)] = (std::make_unique<DiskSurface>(std::string(name), StructureTypeID(sid), StructureID(id), Disk(Vector3(x, y, z), r, Vector3(ux, uy, uz))));
            }
         }

         pos = line.find("Domain");
         if (pos == 0)
         {
            pos = line.find("NULL");      // ignore domain errors
            if (pos > 0 && pos != std::string::npos) continue;
            
            int did, multiplicity, shells;
            sscanf(line.c_str(), "Domain=%d\tMulti=%d\tCount=%d", &did, &multiplicity, &shells);
            size_t off = line.find("Shell");
            for (int i = 0; i < shells; ++i)
            {
               size_t off2 = line.find("Shell", off + 1);
               std::string line_shell = line.substr(off, off2 != std::string::npos ? off2 - off - 1 : off2);
               off = off2;

               int sid; uint code;  double r, x, y, z; char s;
               sscanf(line_shell.c_str(), "Shell=%d\tCode=%d\tRadius=%le\tPosition=%le,%le,%le\t%c", &sid, &code, &r, &x, &y, &z, &s);
               uint clrCode = 0;
               switch (static_cast<Shell::Code>(code))
               {
               case Shell::Code::INIT: clrCode = 0; break;
               case Shell::Code::NORMAL: clrCode = multiplicity; break;
               case Shell::Code::MULTI: clrCode = 3; break;        // multiplicity is not 3, but number of particles.
               default: break;
               }
               if (s == 'S')
                  domains.emplace_back(std::make_pair(Shell(DomainID(did), Sphere(Vector3(x, y, z), r), static_cast<Shell::Code>(code)), clrCode));
               else
               {
                  double ux, uy, uz, hl;
                  sscanf(line_shell.c_str(), "%*s\t%*s\t%*s\t%*s\t%*s\tUnitZ=%le,%le,%le\tHalfLength=%le", &ux, &uy, &uz, &hl);
                  domains.emplace_back(std::make_pair(Shell(DomainID(did), Cylinder(Vector3(x, y, z), r, Vector3(ux, uy, uz), hl), static_cast<Shell::Code>(code)), clrCode));
               }
            }
         }

         pos = line.find("End");
         if (pos == 0)
         {
            particles_.swap(particles);
            domains_.swap(domains);
            structures_.swap(structures);
            time_ = time;
            num_steps_ = num_steps;
            world_size_ = world_size;
            matrix_size_ = matrix_size;

#if defined(_MSC_VER)
            WIN32_FILE_ATTRIBUTE_DATA fileInfoEx;
            GetFileAttributesExA(filename, GET_FILEEX_INFO_LEVELS::GetFileExInfoStandard, &fileInfoEx);
            last_write_hi_ = fileInfoEx.ftLastWriteTime.dwHighDateTime;
            last_write_lo_ = fileInfoEx.ftLastWriteTime.dwLowDateTime;
#else
            // TODO
#endif
            std::strcpy(szFilename, filename);
            sections = section_count;
            break;
         }
      }
      infile.close();
      return true;

   }
   catch (std::out_of_range oor)
   {
      infile.close();
      return false;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

int ExtSim::SelectSection(int selSection)
{
   if (section == selSection) return section;
   if (selSection >= 0 && selSection < static_cast<int>(sections))
   {
      section = selSection;
      readSimFile(szFilename);
      return section;
   }
   return std::min(0, std::max(static_cast<int>(sections), section));
}

// --------------------------------------------------------------------------------------------------------------------------------

bool ExtSim::refresh()
{
#if defined(_MSC_VER)
   WIN32_FILE_ATTRIBUTE_DATA fileInfoEx;
   GetFileAttributesExA(szFilename, GET_FILEEX_INFO_LEVELS::GetFileExInfoStandard, &fileInfoEx);
   if (last_write_hi_ != fileInfoEx.ftLastWriteTime.dwHighDateTime || last_write_lo_ != fileInfoEx.ftLastWriteTime.dwLowDateTime)
      return readSimFile(szFilename);
#else
   // TODO
#endif

   return false;
}

// --------------------------------------------------------------------------------------------------------------------------------
