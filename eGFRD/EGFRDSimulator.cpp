#include "DefsEgfrd.hpp"
#include "EGFRDSimulator.hpp"
#include <fstream>

// --------------------------------------------------------------------------------------------------------------------------------

void format_shape(Structure* st, std::ofstream& f)
{
   auto *box = dynamic_cast<CuboidalRegion*>(st);
   if (box != nullptr) f << "CUBE" << "\tPosition=" << box->position().X() << "," << box->position().Y() << ", " << box->position().Z() << "\tHalfExtent=" << box->shape().half_extent().X() << "," << box->shape().half_extent().Y() << "," << box->shape().half_extent().Z() << "\n";

   auto *plane = dynamic_cast<PlanarSurface*>(st);
   if (plane != nullptr) f << "PLANE" << "\tPosition=" << plane->position().X() << "," << plane->position().Y() << ", " << plane->position().Z() << "\ttHalfExtent=" << plane->shape().half_extent().X() << "," << plane->shape().half_extent().Y() << "\tUnitX=" << plane->shape().unit_x().X() << "," << plane->shape().unit_x().Y() << "," << plane->shape().unit_x().Z() << "\tUnitY=" << plane->shape().unit_y().X() << "," << plane->shape().unit_y().Y() << "," << plane->shape().unit_y().Z() << "\tOneSide=" << plane->shape().is_one_sided() << "\n";

   auto *disk = dynamic_cast<DiskSurface *>(st);
   if (disk != nullptr) f << "DISK" << "\tPosition=" << disk->position().X() << "," << disk->position().Y() << ", " << disk->position().Z() << "\tRadius=" << disk->shape().radius() << "\tUnitZ=" << disk->shape().unit_z().X() << "," << disk->shape().unit_z().Y() << "," << disk->shape().unit_z().Z() << "\n";

   auto *sphere = dynamic_cast<SphericalSurface *>(st);
   if (sphere != nullptr) f << "SPHERE" << "\tPosition=" << sphere->position().X() << "," << sphere->position().Y() << ", " << sphere->position().Z() << "\tRadius=" << sphere->shape().radius() << "\n";

   auto *cyl = dynamic_cast<CylindricalSurface*>(st);
   if (cyl != nullptr) f << "CYL" << "\tPosition=" << cyl->position().X() << "," << cyl->position().Y() << ", " << cyl->position().Z() << "\tRadius=" << cyl->shape().radius() << "\tUnitZ=" << cyl->shape().unit_z().X() << "," << cyl->shape().unit_z().Y() << "," << cyl->shape().unit_z().Z() << "\tHalfLength=" << cyl->shape().half_length() << "\n";
}

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT void EGFRDSimulator::dump(std::string filename, bool append) const
{
   // dump to file in gfrdVisualizer format!

   std::ofstream f;
   f.open(filename, std::ios::out | (append ? std::ofstream::app : std::ofstream::trunc));
   f << std::setprecision(12) << std::scientific;

   f << "sim_step=" << num_steps_ << "\tsim_time=" << time_ << "\t";
   f << "world_size=" << world_.world_size().X() << "\tmatrix_size=" << world_.matrix_size()[0] << "\n";

   for (auto& st : world_.get_structures())
   {
      f << "Structure=" << static_cast<idtype>(st->id()) << "\tName=" << st->name() << "\tSid=" << static_cast<idtype>(st->sid()) << "\tShape=";
      format_shape(st.get(), f);
   }

   for (auto& p : world_.get_particles())
      f << "Particle=" << static_cast<idtype>(p.first) << "\tRadius=" << p.second.radius() << "\tPosition=" << p.second.position().X() << "," << p.second.position().Y() << "," << p.second.position().Z() << "\tSid=" << static_cast<idtype>(p.second.sid()) << "\n";

   for (auto& d : domains_)
   {
      f << "Domain=" << static_cast<idtype>(d.first);
      if (d.second == nullptr) { f << "\tNULL" << "\n"; continue; } // this should not exist but it does :-(
      f << "\tMulti=" << static_cast<uint>(d.second->multiplicity()) << "\tCount=" << d.second->num_shells();
      for (auto &s : d.second->get_shell_list())
      {
         f << "\tShell=" << static_cast<idtype>(s.first) << "\tCode=" << static_cast<uint>(s.second.get().code()) << "\tRadius=";
         switch (s.second.get().shape())
         {
         case Shell::Shape::SPHERE: f << s.second.get().get_sphere().radius() << "\tPosition=" << s.second.get().position().X() << "," << s.second.get().position().Y() << "," << s.second.get().position().Z() << "\tS"; break;
         case Shell::Shape::CYLINDER: f << s.second.get().get_cylinder().radius() << "\tPosition=" << s.second.get().position().X() << "," << s.second.get().position().Y() << "," << s.second.get().position().Z() << "\tC" << "\tUnitZ=" << s.second.get().get_cylinder().unit_z().X() << "," << s.second.get().get_cylinder().unit_z().Y() << "," << s.second.get().get_cylinder().unit_z().Z() << "\tHalfLength=" << s.second.get().get_cylinder().half_length(); break;
         }
      }
      f << "\n";
   }

   f << "End\n";
   f.close();
}

// --------------------------------------------------------------------------------------------------------------------------------
