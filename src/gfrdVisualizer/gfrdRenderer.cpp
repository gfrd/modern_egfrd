#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4505)
#endif
#include <GL/freeglut.h>
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

#include "gfrdRenderer.hpp"
#include "Vector3.hpp"
#include "Matrix4.hpp"
#include "CameraController.hpp"
#include "ExternalSim.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

static Logger& _log = Log("Visualize");

// --------------------------------------------------------------------------------------------------------------------------------

uint colors[] = { 0x000000,
   0x308dd4,0x8dd430,0x66ffff,0x9966ff,
   0x9999ff,0xbfb8af,0xffc89a,0x5bf0a9,
   0xc3eee7,0x6ab3a8,0x5c3170,0xfcac3a,
   0x7eeda2,0x6aec9a,0x008b8b,0xfcf0e5,
   0x397fdb,0xffc89a,0xffe69a,0xff34b3,
   0xfcac3a,0xff4e50,0x4ec5ff,0x6aec9a,
   0x9344bb,0xeed5b7,0xfef1b5,0x990000,
   0xfcac3a,0xdc5757,0x452945,0x412d61,
   0x008b8b,0xdae5e8,0xf6ba00,0xff34b3,
   0x808000,0x838b00,0x008000,
};

// --------------------------------------------------------------------------------------------------------------------------------

GLfloat  cfgStructureLineWidth = 1.2f;
int      cfgStructureCyclinderSides = 32;
int      cfgStructureSphereSlices = 32;
int      cfgStructureSphereStacks = 12;
int      cfgParticleSphereSlices = 12;
int      cfgParticleSphereStacks = 6;
uint     cfgWorldColor = 0x0336699;

uint     cfgShellColor = 0x0F00000;
uint     cfgSelectedShellColor = 0x0F0E0D0;
uint     cfgSelectedStructureColor = 0x0D0E0F0;
uint     cfgSelectedParticleColor = 0x0D4AF37;

uint     cfgShellColors[] = { 0x0A050F0 , 0x0F0F000 , 0x000F000 , 0x0F00000, 0x0D4AF37 };    // INIT(violet), SINGLE(yellow), PAIR(green), MULTI(red), SELECTED(gold)

// --------------------------------------------------------------------------------------------------------------------------------

void RenderGFRD(const EGFRDSimulator& s, const CameraController& cam, int showpid, bool drawShells, StructureTypeID selStructure, ParticleID selParticle, DomainID selDomain, std::array<uint, 3>& domain_type_count)
{
   const World& w = s.world();
   drawMatrix(w.cell_size(), w.matrix_size());

   // particles
   auto v = cam.calculate(Vector3::uy).normal();
   glEnable(GL_LIGHTING);
   for (auto p : w.get_particles())
      drawParticle(p, selParticle, showpid, v, &s.world());

   drawAllStructures(w.get_structures(), cam, selStructure);

   if (!drawShells) return;

   std::vector<std::pair<std::reference_wrapper<const Shell>, uint>> shells;
   for (auto& domain : s.get_domains())
   {
      domain_type_count[static_cast<int>(domain.second->multiplicity()) - 1]++;
      if (domain.second->num_shells() == 1)
      {
         uint clrCode = selDomain == domain.first ? 4 : domain.second->get_shell().second.get().code() == Shell::Code::INIT ? 0 : static_cast<int>(domain.second->multiplicity());
         shells.emplace_back(std::make_pair(domain.second->get_shell().second, clrCode));
      }
      else
      {
         uint clrCode = selDomain == domain.first ? 4u : 3u;
         for (const auto& sis : domain.second->get_shell_list())
            shells.emplace_back(std::make_pair(sis.second, clrCode));
      }
   }

   const auto u = cam.calculate(Vector3::ux).normal();
   drawAllShells(shells, u);
}

// --------------------------------------------------------------------------------------------------------------------------------

void RenderGFRDExtern(const ExtSim& w, const CameraController& cam, int showpid, bool drawShells, uint selStructure, uint selParticle, std::array<uint, 3>& domain_type_count)
{
   drawMatrix(w.cell_size(), w.matrix_size());
   drawBox(Box(w.world_size() / 2, w.world_size() / 2));

   // particles
   auto v = cam.calculate(Vector3::uy).normal();
   glEnable(GL_LIGHTING);
   for (auto p : w.get_particles())
      drawParticle(p, ParticleID(selParticle), showpid, v, nullptr);

   drawAllStructures(w.get_structures(), cam, StructureTypeID(selStructure));

   if (!drawShells) return;

   std::vector<std::pair<std::reference_wrapper<const Shell>, uint>> shells;
   DomainID did;
   for (const auto& sp : w.get_domains())
   {
      if (did != sp.first.did()) domain_type_count[sp.second > 0 ? sp.second - 1 : 0]++;    // count init as single, don't overcount multi's
      did = sp.first.did();
      shells.emplace_back(std::make_pair(std::cref(sp.first), sp.second));
   }

   auto u = cam.calculate(Vector3::ux).normal();
   drawAllShells(shells, u);
}

// --------------------------------------------------------------------------------------------------------------------------------

void drawAllStructures(StructureContainer::structures_range s, const CameraController& cam, StructureTypeID selStructure)
{
   UNUSED(cam);
   for (auto st : s)
   {
      auto *box = dynamic_cast<CuboidalRegion*>(st.get());
      if (box != nullptr) drawBox(box->shape());

      auto *plane = dynamic_cast<PlanarSurface*>(st.get());
      if (plane != nullptr) drawPlane(plane->shape(), plane->sid(), plane->sid() == selStructure);

      auto *disk = dynamic_cast<DiskSurface *>(st.get());
      if (disk != nullptr) drawDisk(disk->shape(), disk->sid(), disk->sid() == selStructure);

      auto *sphere = dynamic_cast<SphericalSurface *>(st.get());
      if (sphere != nullptr) drawSphere(sphere->shape(), sphere->sid(), sphere->sid() == selStructure);

      auto *cyl = dynamic_cast<CylindricalSurface*>(st.get());
      if (cyl != nullptr) drawCylinder(cyl->shape(), cyl->sid(), cyl->sid() == selStructure);
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

void drawAllShells(const std::vector<std::pair<std::reference_wrapper<const Shell>, uint>>& shells, const Vector3& v)
{
   for (const auto& spair : shells)
   {
      const auto& shell = spair.first.get();

      auto color = cfgShellColors[spair.second];
      glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);

      glPushMatrix();
      glTranslated(shell.position().X(), shell.position().Y(), shell.position().Z());
      glEnable(GL_LIGHTING);

      if (shell.shape() == Shell::Shape::SPHERE)
      {
         glRotated(90.0, 1.0, 0, 0);      // rotate around X draw Sphere upward (nicer)
         glutWireSphere(shell.get_sphere().radius(), 24, 8);
      }
      if (shell.shape() == Shell::Shape::CYLINDER)
      {
         auto c = shell.get_cylinder();
         auto m = Matrix4::createRotationAB(Vector3::uz, c.unit_z());
         double mm[16];
         glMultMatrixd(m.getGLMatrix(mm));
         glTranslated(0, 0, -c.half_length());
         glutWireCylinder(c.radius(), 2 * c.half_length(), 24, 1);
      }

      glDisable(GL_LIGHTING);
      glPopMatrix();

      glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);
      auto x1 = shell.position() + 1.05 * (shell.shape() == Shell::Shape::SPHERE ? shell.get_sphere().radius() : shell.get_cylinder().radius()) * v;
      glRasterPos3d(x1.X(), x1.Y(), x1.Z());
      glutBitmapString(GLUT_BITMAP_8_BY_13, reinterpret_cast<const unsigned char*>(static_cast<std::string>(make_string() << static_cast<idtype>(shell.did())).c_str()));
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

void drawMatrix(double cs, std::array<uint, 3> ms)
{
   uint color = cfgWorldColor;
   glLineWidth(1.0);
   glLineStipple(1, 0x3333);
   glEnable(GL_LINE_STIPPLE);
   glDisable(GL_LIGHTING);
   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);
   glBegin(GL_LINES);
   for (uint i = 1; i < ms[0]; ++i)
   {
      glVertex3d(i*cs, 0, 0);
      glVertex3d(i*cs, 0, cs*ms[2]);
   }
   glBegin(GL_LINES);
   for (uint i = 1; i < ms[2]; ++i)
   {
      glVertex3d(0, 0, i*cs);
      glVertex3d(cs*ms[0], 0, i*cs);
   }
   glEnd();
   glEnable(GL_LIGHTING);
   glDisable(GL_LINE_STIPPLE);
}

// --------------------------------------------------------------------------------------------------------------------------------

void drawParticle(const Particle::particle_id_pair& pip, ParticleID selParticle, int showID, const Vector3& v, const World* world)
{
   glEnable(GL_LIGHTING);

   glPushMatrix();
   glTranslated(pip.second.position().X(), pip.second.position().Y(), pip.second.position().Z());

   auto color = selParticle == pip.first ? cfgSelectedParticleColor : colors[static_cast<idtype>(pip.second.sid())];
   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);
   glutSolidSphere(pip.second.radius(), cfgParticleSphereSlices, cfgParticleSphereStacks);
   glPopMatrix();

   if (!showID) return;

   auto x1 = pip.second.position() + 1.05 * pip.second.radius() * v;
   glDisable(GL_LIGHTING);
   glRasterPos3d(x1.X(), x1.Y(), x1.Z());

   if (showID == 1)
      glutBitmapString(GLUT_BITMAP_8_BY_13, reinterpret_cast<const unsigned char*>(static_cast<std::string>(make_string() << static_cast<idtype>(pip.first)).c_str()));
   else
   {
      auto &name = world ? world->get_species(pip.second.sid()).name() : make_string() << static_cast<idtype>(pip.second.sid());
      glutBitmapString(GLUT_BITMAP_8_BY_13, reinterpret_cast<const unsigned char*>(name.c_str()));
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

void drawBox(const Box& b)
{
   glDisable(GL_LIGHTING);
   uint color = cfgWorldColor;
   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);

   glPushMatrix();
   glTranslated(b.half_extent().X(), b.half_extent().Y(), b.half_extent().Z());
   glScaled(b.half_extent().X(), b.half_extent().Y(), b.half_extent().Z());

   glLineWidth(1.0);
   glLineStipple(1, 0x3333);
   glEnable(GL_LINE_STIPPLE);
   glutWireCube(2.0);
   glDisable(GL_LINE_STIPPLE);
   glPopMatrix();
}

// --------------------------------------------------------------------------------------------------------------------------------

void drawPlane(const Plane& p, const StructureTypeID& sid, bool select)
{
   auto color = select ? cfgSelectedStructureColor : colors[static_cast<idtype>(sid)];
   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);

   glPushMatrix();
   glTranslated(p.position().X(), p.position().Y(), p.position().Z());

   Vector3 lx = p.half_extent().X() * p.unit_x();
   Vector3 ly = p.half_extent().Y() * p.unit_y();

   glLineWidth(cfgStructureLineWidth);
   glDisable(GL_LIGHTING);
   glDisable(GL_CULL_FACE);
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
   glBegin(GL_LINE_STRIP);
   glVertex3(lx + ly);
   glVertex3(lx - ly);
   glVertex3(-lx - ly);
   glVertex3(-lx + ly);
   glVertex3(lx + ly);
   glEnd();

   if (!select)
   {
      glEnable(GL_LIGHTING);
      if (p.is_one_sided()) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glScalef(0.998f, 0.998f, 0.998f); // slightly smaller, don't overdraw the lines
      glBegin(GL_TRIANGLE_STRIP);
      glVertex3(lx + ly);
      glVertex3(-lx + ly);
      glVertex3(lx - ly);
      glVertex3(-lx - ly);
      glEnd();
   }

   glEnable(GL_CULL_FACE);
   glPopMatrix();
}

void drawDisk(const Disk& d, const StructureTypeID& sid, bool select)
{
   auto color = select ? cfgSelectedStructureColor : colors[static_cast<idtype>(sid)];
   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);

   glPushMatrix();
   glTranslated(d.position().X(), d.position().Y(), d.position().Z());
   auto m = Matrix4::createRotationAB(Vector3::uz, d.unit_z());
   double mm[16];
   glMultMatrixd(m.getGLMatrix(mm));

   glDisable(GL_LIGHTING);
   glLineWidth(cfgStructureLineWidth);
   glutWireCylinder(d.radius(), 0.0, cfgStructureCyclinderSides, 1);
   if (!select)
   {
      glEnable(GL_LIGHTING);
      glutSolidCylinder(d.radius(), 0.0, cfgStructureCyclinderSides, 1);
   }
   glPopMatrix();
}

void drawSphere(const Sphere& s, const StructureTypeID& sid, bool select)
{
   auto color = select ? cfgSelectedStructureColor : colors[static_cast<idtype>(sid)];
   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);

   glPushMatrix();
   glTranslated(s.position().X(), s.position().Y(), s.position().Z());
   glRotated(90.0, 1.0, 0, 0);      // rotate around X draw Sphere upward (nicer)

   glDisable(GL_LIGHTING);
   glLineWidth(cfgStructureLineWidth);
   glutWireSphere(s.radius(), cfgStructureSphereSlices, cfgStructureSphereStacks);
   if (!select)
   {
      glEnable(GL_LIGHTING);
      glutSolidSphere(s.radius(), cfgStructureSphereSlices, cfgStructureSphereStacks);
   }
   glPopMatrix();
}

void glxWireCylinder(double radius, double length, int sides)
{
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
   glBegin(GL_LINES);
   for (int i = 0; i <= sides; i++)
   {
      double alpha = M_PI * 2.0 * i / sides;
      glVertex3d(radius*std::sin(alpha), radius*std::cos(alpha), -length);
      glVertex3d(radius*std::sin(alpha), radius*std::cos(alpha), length);
   }
   glEnd();
}

void glxSolidCylinder(double radius, double length, int sides)
{
   glDisable(GL_CULL_FACE);
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   glBegin(GL_TRIANGLE_STRIP);
   for (int i = 0; i <= sides; i++)
   {
      double alpha = M_PI * 2.0 * i / sides;
      glVertex3d(radius*std::sin(alpha), radius*std::cos(alpha), -length);
      glVertex3d(radius*std::sin(alpha), radius*std::cos(alpha), length);
   }
   glEnd();
   glEnable(GL_CULL_FACE);
}

void drawCylinder(const Cylinder& c, const StructureTypeID& sid, bool select)
{
   auto color = select ? cfgSelectedStructureColor : colors[static_cast<idtype>(sid)];
   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);

   glPushMatrix();
   glTranslated(c.position().X(), c.position().Y(), c.position().Z());
   auto m = Matrix4::createRotationAB(Vector3::uz, c.unit_z());
   double mm[16];
   glMultMatrixd(m.getGLMatrix(mm));

   glDisable(GL_LIGHTING);
   glLineWidth(cfgStructureLineWidth);
   glxWireCylinder(c.radius(), c.half_length(), cfgStructureCyclinderSides);
   if (!select)
   {
      glEnable(GL_LIGHTING);
      glxSolidCylinder(c.radius(), c.half_length(), cfgStructureCyclinderSides);
   }
   glPopMatrix();
}

// --------------------------------------------------------------------------------------------------------------------------------

