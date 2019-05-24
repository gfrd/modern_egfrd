#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4505)
#endif
#include <GL/freeglut.h>
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

#include <Logger.hpp>
#include <EGFRDSimulator.hpp>
#include <Matrix4.hpp>
#include "CameraController.hpp"
#include "gfrdRenderer.hpp"
#include "ScreenShot.hpp"
#include "ExternalSim.hpp"
#include "Persistence.hpp"
#include "../Common/getoptions.hpp"
#include "SimulatorSettings.hpp"

#if defined(_MSC_VER)
#include "resource.h"
#endif

// --------------------------------------------------------------------------------------------------------------------------------

static float ar;
static CameraControllerAnimated cam;
static EGFRDSimulator* sptr;
bool drawShells = false;
bool drawOrig = true;
int drawPid = 0;
bool autoSim = false;
bool autoCheck = false;
bool showHelp = false;
idtype selParticle = 0;
idtype selDomain = 0;
idtype selStructure = 0;
uint selSection = 0;
std::string modelfile;

// --------------------------------------------------------------------------------------------------------------------------------

std::vector<std::pair<bool, std::string>> help = 
{
   std::make_pair(true, "ESC = exit"),
   std::make_pair(true, "H = show help") ,
   std::make_pair(true, "F = full screen") ,
   std::make_pair(true, "S = show shells") ,
   std::make_pair(true, "D = demo rotate") ,
   std::make_pair(true, "P = screen-shot") ,
   std::make_pair(true, "O = draw origin") ,
   std::make_pair(true, "I = particle ID/SID") ,
   std::make_pair(false, ". = sim step") ,
   std::make_pair(false, "/ = check sim") ,
   std::make_pair(false, "A = auto") ,
   std::make_pair(false, "B = burst all") ,
   std::make_pair(true, "</> = select structure") ,
   std::make_pair(true, "+/- = select particle") ,
   std::make_pair(true, "g = goto particle") ,
   std::make_pair(true, "0 = deselect all") ,
   std::make_pair(true, "[/] = select domain") ,
   std::make_pair(true, "j = goto domain") ,
   std::make_pair(true, "C = center world"),
};

// --------------------------------------------------------------------------------------------------------------------------------

size_t get_help_width()
{
   size_t width = 0;
   for (auto& txt : help)
   {
      size_t size = glutBitmapWidth(GLUT_BITMAP_9_BY_15, 'c') * txt.second.length();
      if (txt.first || sptr != nullptr) width = std::max(width, size);
   }
   return width;
}

// --------------------------------------------------------------------------------------------------------------------------------

void check_sim()
{
   try
   {
      sptr->check();
   }
   catch (const gfrd_exception& ex)
   {
      Log("Visualize").fatal() << "Check failed: time: " << sptr->time() << ", step: " << sptr->num_steps() << ", fault: " << ex.what();
      autoSim = false;
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

void handleKeyboard(unsigned char cChar, int nMouseX, int nMouseY)
{
   UNUSED(nMouseX, nMouseY);
   try
   {
      switch (cChar)
      {
         case 27: glutLeaveMainLoop(); break;
         case 'f': case 'F': glutFullScreenToggle(); break;
         case 'h': case 'H': showHelp = !showHelp; glutPostRedisplay(); break;
         case 's': case 'S': drawShells = !drawShells; glutPostRedisplay(); break;
         case 'd': case 'D': cam.set_demo(!cam.demo()); break;
         case 'p': case 'P': screenshot(static_cast<int>(sptr ? sptr->num_steps() : extSim.num_steps())); break;
         case 'o': case 'O': drawOrig = !drawOrig; glutPostRedisplay(); break;
         case 'i': case 'I': drawPid++; if (drawPid == 3) drawPid = 0; glutPostRedisplay(); break;

         case '.': if (sptr) sptr->step(); glutPostRedisplay(); break;
         case '/': if (sptr) check_sim(); break;
         case 'a': case 'A': autoSim = !autoSim; autoCheck = (cChar == 'a'); glutPostRedisplay(); break;
         case 'b': case 'B': if (sptr) sptr->burst_all(); glutPostRedisplay(); break;

         case '>': selStructure++; glutPostRedisplay(); break;
         case '<': if (selStructure > 0) selStructure--; glutPostRedisplay(); break;

         case '+': selParticle++; glutPostRedisplay(); break;
         case '-': if (selParticle > 0) { selParticle--; glutPostRedisplay(); } break;
         case 'g': case 'G': std::cout << "Enter ParticleID: ";  std::cin >> selParticle; if (selParticle > 0) { glutPostRedisplay(); } break;

         case ']': selDomain++; glutPostRedisplay(); break;
         case '[': if (selDomain > 0) { selDomain--; glutPostRedisplay(); } break;
         case 'j': case 'J': std::cout << "Enter DomainID: ";  std::cin >> selDomain; if (selDomain > 0) { glutPostRedisplay(); } break;

         case '0': selParticle = 0; selDomain = 0;  glutPostRedisplay(); break;
         case 'c': case 'C': selParticle = 0;
            if (sptr)
            {
               const auto& w = sptr->world();
               double side = 0.5 / std::max(w.world_size().X(), std::max(w.world_size().Y(), w.world_size().Z()));      // scale largest side to GLunits 1.0
               cam.lookTo(side * w.world_size());
            }
            else
            {
               double side = 0.5 / std::max(extSim.world_size().X(), std::max(extSim.world_size().Y(), extSim.world_size().Z()));      // scale largest side to GLunits 1.0
               cam.lookAt(side * extSim.world_size());
            }
            glutPostRedisplay();
            break;

         case 'q': case 'Q': selSection = extSim.SelectSection(--selSection); glutPostRedisplay(); break;
         case 'w': case 'W': selSection = extSim.SelectSection(++selSection); glutPostRedisplay(); break;

            //case 'n': case 'j': case 'm': case 'k': case 'l': case ',':    // for overlap checking , move particle1 around!
            //{
            //   const double step = 0.05;
            //   auto& w = sptr->world();
            //   auto pip = w.get_particle(ParticleID(1));
            //   if (cChar == 'n') pip.second.position() += Vector3(-step, 0, 0);
            //   if (cChar == 'j') pip.second.position() += Vector3(+step, 0, 0);
            //   if (cChar == 'm') pip.second.position() += Vector3(0, -step, 0);
            //   if (cChar == 'k') pip.second.position() += Vector3(0, +step, 0);
            //   if (cChar == 'l') pip.second.position() += Vector3(0, 0, -step);
            //   if (cChar == ',') pip.second.position() += Vector3(0, 0, +step);
            //   const_cast<World&>(w).update_particle(pip);
            //   auto s = w.get_structure(StructureID(2));
            //   auto cyl = dynamic_cast<CylindricalSurface*>(s.get());
            //   const Cylinder shape = cyl->shape();
            //   auto ovl = w.check_particle_overlap(shape);
            //   if (ovl.size() > 0)
            //      for (auto& i : ovl)
            //         std::cout << i.first.first << " overlaps " << i.second << " : " << (-i.second > 2 * pip.second.radius() ? "IN" : "OUT") << std::endl;
            //   else
            //      std::cout << "no overlaps" << std::endl;
            //   glutPostRedisplay();
            //} break;

         case '1': cam.set_angles(0.0, M_PI / 2); glutPostRedisplay(); break;      // TOP
         case '2': cam.set_angles(0.0, -M_PI / 2); glutPostRedisplay(); break;     // BOTTOM
         case '3': cam.set_angles(0.0, 0.0); glutPostRedisplay(); break;         // FRONT
         case '4': cam.set_angles(M_PI / 2, 0.0); glutPostRedisplay(); break;      // LEFT
         case '5': cam.set_angles(M_PI, 0.0); glutPostRedisplay(); break;        // BACK
         case '6': cam.set_angles(-M_PI / 2, 0.0); glutPostRedisplay(); break;     // RIGHT

         default: break;
      }
   }
   catch (const std::exception&)
   {
      abort();
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

void handleReshape(int width, int height)
{
   ar = static_cast<float>(width) / static_cast<float>(height);
   glViewport(0, 0, width, height);
   cam.handleReshape(width, height);
}

// --------------------------------------------------------------------------------------------------------------------------------

void handleIdle()
{
   try
   {
      cam.handleIdle();
      if (autoSim && sptr)
      {
         if (!sptr->step()) autoSim = false;
         if (autoCheck) check_sim();
         glutPostRedisplay();
      }

      if (extSim.active() && extSim.refresh()) glutPostRedisplay();
   }
   catch (const std::exception& ex)
   {
      printf("Unhandled exception: %s", ex.what());
      abort();
   }
}

// --------------------------------------------------------------------------------------------------------------------------------

void drawPoint(Vector3 point, double size, uint color)
{
   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);

   glLineWidth(2.f);
   glDisable(GL_LIGHTING);

   glBegin(GL_LINES);
   glVertex3(point - Vector3(size, 0, 0));
   glVertex3(point + Vector3(size, 0, 0));
   glVertex3(point - Vector3(0, size, 0));
   glVertex3(point + Vector3(0, size, 0));
   glVertex3(point - Vector3(0, 0, size));
   glVertex3(point + Vector3(0, 0, size));
   glEnd();
}

void drawArrow(Vector3 p1, Vector3 p2, uint color, double size)
{
   auto l = (p2 - p1).length();
   auto d = (p2 - p1).normal();
   auto t = p2 - Vector3::multiply(d, size);

   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();

   glTranslated(t.X(), t.Y(), t.Z());
   auto m = Matrix4::createRotationAB(Vector3::uz, d);
   double mm[16];
   glMultMatrixd(m.getGLMatrix(mm));

   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);
   glBegin(GL_LINES);
   glVertex3d(0, 0, -(l - size));
   glVertex3d(0, 0, 0);
   glEnd();

   glutWireCone(size / 8.f, size, 6, 1);
   glPopMatrix();
}


void drawOrigin()
{
   const uint color = 0x00F0EEDF;
   auto scale = std::abs(cam.distance() / 20.f);
   auto pos = Vector3(-0.2 * scale, -0.2 * scale, -0.2 * scale);
   drawArrow(pos, pos + Vector3::multiply(Vector3::ux, scale), 0x00FF0000, scale / 8);
   drawArrow(pos, pos + Vector3::multiply(Vector3::uy, scale), 0x0000FF00, scale / 8);
   drawArrow(pos, pos + Vector3::multiply(Vector3::uz, scale), 0x000000FF, scale / 8);

   glColor3ub((color >> 16) & 0xFF, (color >> 8) & 0xFF, (color >> 0) & 0xFF);
   auto x1 = pos + Vector3::multiply(Vector3::ux, 1.1* scale);
   glRasterPos3d(x1.X(), x1.Y(), x1.Z());
   glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>("X"));
   auto y1 = pos + Vector3::multiply(Vector3::uy, 1.1* scale);
   glRasterPos3d(y1.X(), y1.Y(), y1.Z());
   glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>("Y"));
   auto z1 = pos + Vector3::multiply(Vector3::uz, 1.1* scale);
   glRasterPos3d(z1.X(), z1.Y(), z1.Z());
   glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>("Z"));

   if (cam.drag())
      drawPoint(cam.center(), scale / 10, 0x00A0EE2F);
}

// --------------------------------------------------------------------------------------------------------------------------------

void handleDisplay()
{
   try
   {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluPerspective(25.0f, static_cast<GLfloat>(ar), 0.01f, 200.0f);

      Vector3 ws = sptr ? sptr->world().world_size() : extSim.world_size();

      double side = 1.0 / std::max(ws.X(), std::max(ws.Y(), ws.Z()));      // scale largest side to GLunits 1.0
      Matrix4 scale = Matrix4::createScale(side * Vector3::one);
      double mscale[16];
      scale.getGLMatrix(mscale);

      // follow particle
      bool selParExists = false;
      bool selDomExists = false;

      if (selParticle > 0)
      {
         auto selPid = ParticleID(selParticle);
         if (sptr != nullptr && sptr->world().has_particle(selPid))
         {
            auto p = sptr->world().get_particle(selPid);
            cam.lookTo(side * p.second.position());
            selParExists = true;
         }
         else
         {
            const auto& particles = extSim.get_particles();
            for (const auto& p : particles)
            {
               if (p.first == selPid)
               {
                  cam.lookTo(side * p.second.position());
                  selParExists = true;
                  break;
               }
            }
         }
      }
      else
         if (selDomain > 0)
         {
            auto selD = DomainID(selDomain);
            if (sptr != nullptr && sptr->has_domain(selD))
            {
               auto& d = sptr->get_domain(selD);
               Vector3 look(0, 0, 0);
               if (d.multiplicity() == Domain::Multiplicity::MULTI)
               {
                  for (auto &s : d.get_shell_list())
                     look += s.second.get().position();
                  look /= static_cast<double>(d.num_shells());
               }
               else
                  look = d.get_shell().second.get().position();
               cam.lookTo(side * look);
               selDomExists = true;
            }
         }

      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      cam.look();
      if (drawOrig) drawOrigin();

      glMultMatrixd(mscale);

      std::array<uint, 3> domain_type_count{ 0,0,0 };
      if (sptr != nullptr)
         RenderGFRD(*sptr, cam, drawPid, drawShells, StructureTypeID(selStructure), ParticleID(selParticle), DomainID(selDomain), domain_type_count);
      else
         RenderGFRDExtern(extSim, cam, drawPid, drawShells, selStructure, selParticle, domain_type_count);

      {
         glMatrixMode(GL_PROJECTION);
         glLoadIdentity();
         glMatrixMode(GL_MODELVIEW);
         glLoadIdentity();

         int viewport[4];
         glGetIntegerv(GL_VIEWPORT, viewport);
         glOrtho(0, viewport[2], 0, viewport[3], -1, 1);

         glDisable(GL_LIGHTING);
         glColor3d(.8, .9, .7);
         int y = 3;
         glRasterPos2i(2, y - glutBitmapHeight(GLUT_BITMAP_9_BY_15) + viewport[3]);
         std::string msg1 = make_string() << "Time: " << (sptr ? sptr->time() : extSim.time()) << " Step:" << (sptr ? sptr->num_steps() : extSim.num_steps());
         glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>(msg1.c_str()));

         if (drawShells)
         {
            y -= 18;
            glRasterPos2i(2, y - glutBitmapHeight(GLUT_BITMAP_9_BY_15) + viewport[3]);
            std::string msg2 = make_string() << "S: " << domain_type_count[0] << " P: " << domain_type_count[1] << " M: " << domain_type_count[2];
            glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>(msg2.c_str()));
         }

         if (selParticle > 0)
         {
            y -= 18;
            glRasterPos2i(2, y - glutBitmapHeight(GLUT_BITMAP_9_BY_15) + viewport[3]);
            std::string msg2 = make_string() << "sel PID: " << selParticle << (selParExists ? "." : "!");
            glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>(msg2.c_str()));
         }

         if (selDomain > 0)
         {
            y -= 18;
            glRasterPos2i(2, y - glutBitmapHeight(GLUT_BITMAP_9_BY_15) + viewport[3]);
            std::string msg2 = make_string() << "sel DID: " << selDomain << (selDomExists ? "." : "!");
            glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>(msg2.c_str()));
         }

         if (selStructure > 0)
         {
            y -= 18;
            glRasterPos2i(2, y - glutBitmapHeight(GLUT_BITMAP_9_BY_15) + viewport[3]);
            std::string msg2 = make_string() << "sel SID: " << selStructure;
            glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>(msg2.c_str()));
         }

         if (showHelp)
         {
            y = 3;
            static size_t help_width = 0;
            if (help_width == 0) help_width = get_help_width();
            for (auto& txt : help)
            {
               if (txt.first || sptr != nullptr) 
               {
                  glRasterPos2i(viewport[2] - static_cast<GLint>(help_width), y - glutBitmapHeight(GLUT_BITMAP_9_BY_15) + viewport[3]);
                  glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>(txt.second.c_str()));
                  y -= 18;
               }
            }
         }

         if (sptr)
         {
            glRasterPos2i(2, 2);    // Copy Numbers
            std::stringstream line;
            for (auto& s : sptr->world().get_species())
               line << s.second.name() << ": " << sptr->world().get_particle_ids(s.first).size() << " ";
            glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>(line.str().c_str()));
         }
         else
         {
            glRasterPos2i(2, 2);    // ext sim sections
            std::stringstream line;
            line << extSim.Section() << "/" << extSim.Sections();
            glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>(line.str().c_str()));
         }

      }

      glutSwapBuffers();
   }
   catch (const std::exception&)
   {
      abort();
   }

   // Record all frames for movie
   //screenshot(number++);
}

// --------------------------------------------------------------------------------------------------------------------------------

void handleMouse(int button, int updown, int x, int y)
{
   cam.handleMouse(button, updown, x, y);
}

// --------------------------------------------------------------------------------------------------------------------------------

void handlePassiveMouseMotion(int x, int y)
{
   cam.handlePasiveMouseMotion(x, y);
}

// --------------------------------------------------------------------------------------------------------------------------------

void handleMouseMotion(int x, int y)
{
   cam.handleMouseMotion(x, y);
}

// --------------------------------------------------------------------------------------------------------------------------------

void handleMouseWheel(int wheel_number, int direction, int x, int y)
{
   cam.handleMouseWheel(wheel_number, direction, x, y);
}

// --------------------------------------------------------------------------------------------------------------------------------

void SetupModel(const std::string& filename, const std::vector < std::string>& arguments, World& world, RandomNumberGenerator& rng, ReactionRuleCollection& rules)
{
   SimulatorSettings settings;

   for (auto& arg : arguments)
      settings.add_variable(arg);

   std::cout << std::setw(14) << "model file = " << filename << "\n";

   std::ifstream stream(filename);
   if (!stream.is_open()) { Log("RunGfrd").fatal() << "Failed to load input file!"; return; }
   stream >> settings;
   stream.close();

   {
      auto &section = settings.getSimulatorSection();
      if (section.seed()) rng.seed(section.seed());
   }

   double world_size;
   {
      auto& section(settings.getWorldSection());
      world_size = section.world_size();
      THROW_UNLESS_MSG(illegal_size, world_size > 0, "invalid world-size!");
      THROW_UNLESS_MSG(illegal_size, world.matrix_size()[0] == section.matrix_space(), "invalid matrix space x-size!");
      THROW_UNLESS_MSG(illegal_size, world.matrix_size()[1] == section.matrix_space(), "invalid matrix space y-size!");
      THROW_UNLESS_MSG(illegal_size, world.matrix_size()[2] == section.matrix_space(), "invalid matrix space z-size!");
   }

   auto& vars = settings.getVariablesSection();

   Model model;
   {
      auto& sections = settings.getSpeciesTypeSections();
      for (auto& section : sections)
         section.create_species(model, vars);
   }

   world.initialize(world_size, model);

   for (auto& section : settings.getReactionRuleSections())
      section.create_reaction_rule(model, rules);

   for (auto& section : settings.getParticlesSections())
      section.add_particles_to_world(model, world, rng);


   // Print !

   std::cout << "\n";

   settings.getVariablesSection().PrintSettings();

   //Simulation::PrintSettings();

   std::cout << "\n";

   for (auto& species : settings.getSpeciesTypeSections())
      species.PrintSettings();

   std::cout << "\n";

   for (auto& rule : settings.getReactionRuleSections())
      rule.PrintSettings();

   std::cout << "\n";

   for (auto& section : settings.getParticlesSections())
      section.PrintSettings();

   std::cout << "\n";

   const auto& cns = settings.getCopyNumbersSection();
   if (cns != nullptr) std::cout << "CopyNumbersSection ignored in visualizer.\n";
   const auto& pps = settings.getParticlePositionsSection();
   if (pps != nullptr) std::cout << "ParticlePositionsSection ignored in visualizer.\n";
   const auto& rrs = settings.getReactionRecordSection();
   if (rrs != nullptr) std::cout << "ReactionRecordSection ignored in visualizer.\n";
   const auto& ps = settings.getProgressSection();
   if (ps != nullptr) std::cout << "ProgressSection ignored in visualizer.\n";

}

// --------------------------------------------------------------------------------------------------------------------------------

void print_usage()
{
   std::cout << "  1) gfrdVisualizer <path-to-model> [ -d var=value ]" << std::endl;
   std::cout << "  2) gfrdVisualizer -r,--resume <path-to-simstate> [glut-options]" << std::endl;
   std::cout << "  3) gfrdVisualizer -c,--crash <path-to-dump> [glut-options]" << std::endl;
   std::cout << "  4) gfrdVisualizer -x,--demo [ N ] [glut-options]" << std::endl;
   std::cout << "        [-h,-?,--help]        Print command line usage information" << std::endl;
   std::cout << "        [-v,--version]        Print version/build information" << std::endl << std::endl;
   std::cout << "1) Start simulation described in model-file." << std::endl;
   std::cout << "2) Load simulator maintenance/state file." << std::endl;
   std::cout << "3) Load simulator crash/dump file." << std::endl;
   std::cout << "4) Start simulation build-in demo (N = 1..3)" << std::endl;
   std::cout << std::endl << std::endl;
}

// --------------------------------------------------------------------------------------------------------------------------------

void initRenderState()
{
   glClearColor(0, 0, 0, 0);
   glEnable(GL_CULL_FACE);
   glCullFace(GL_BACK);

   glEnable(GL_DEPTH_TEST);                        // Enables Depth Testing
   glDepthFunc(GL_LEQUAL);                         // The Type Of Depth Testing To Do
   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);          // Really Nice Perspective Calculations
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   glEnable(GL_LIGHT0);
   glEnable(GL_NORMALIZE);
   glEnable(GL_COLOR_MATERIAL);

   const GLfloat light_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
   const GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
   const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
   const GLfloat light_position[] = { 2.0f, 5.0f, 5.0f, 0.0f };

   glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
   glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

   const GLfloat mat_ambient[] = { 0.7f, 0.7f, 0.7f, 1.0f };
   const GLfloat mat_diffuse[] = { 0.8f, 0.8f, 0.8f, 1.0f };
   const GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
   const GLfloat high_shininess[] = { 20.0f };

   glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
   glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, high_shininess);
}

// --------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
   // local to main, need to stay in scope
   RandomNumberGenerator rng;
   ReactionRuleCollection rules;
   World world;
   EGFRDSimulator s(world, rules, rng);
   sptr = &s;
   std::vector<std::string> model_arguments;

   int localDemoInit = 0;
   try
   {
      int arg_err = -1;
      getoptions args(argc, argv);
      for (size_t i = 0; i < args.size() && arg_err == -1; ++i)
      {
         if (localDemoInit == 0 && args.isvalue_F(i))
         {
            modelfile = args.option(i);
            localDemoInit = 1;
            continue;
         }

         if ((localDemoInit == 0 || localDemoInit == 1) && args.isparam(i) && args.option(i) == "d" && args.isvalue(i + 1))
         {
            model_arguments.emplace_back(args.option(++i));
            continue;
         }

         if (localDemoInit == 0 && (args.option(i) == "r" || args.option(i) == "-resume") && args.isvalue_F(i + 1))
         {
            modelfile = args.option(++i);
            localDemoInit = 2;
            continue;
         }

         if (localDemoInit == 0 && (args.option(i) == "c" || args.option(i) == "-crash") && args.isvalue_F(i + 1))
         {
            modelfile = args.option(++i);
            localDemoInit = 3;
            continue;
         }

         if (localDemoInit == 0 && (args.option(i) == "x" || args.option(i) == "-demo"))
         {
            if (args.isvalue_NP(i + 1)) localDemoInit = 100 + std::stoi(args.option(++i));
            else localDemoInit = 101;
            continue;
         }

         if (args.option(i) == "h" || args.option(i) == "?" || args.option(i) == "-help")
         {
            gfrd_print_header();
            std::cout << "Usage:" << std::endl << std::endl;
            print_usage();
            return 1;
         }

         if (args.option(i) == "v" || args.option(i) == "-version")
         {
            gfrd_print_version();
            return 1;
         }

         // glut-options, ignore those here
         if ((args.option(i) == "display" || args.option(i) == "geometry") && !args.isparam(i + 1)) { i++; continue; }
         if (args.option(i) == "iconic" || args.option(i) == "indirect") continue;
         if (args.option(i) == "gldebug" || args.option(i) == "direct" || args.option(i) == "sync") continue;

         // all that remains is an error
         arg_err = static_cast<int>(i);
      }

      if (arg_err != -1)
      {
         gfrd_print_header();
         std::cout << "ERROR: Unknown or invalid argument: " << (args.isparam(arg_err) ? "-" : "") << args.option(arg_err) << std::endl;
         std::cout << "use --help argument to print usage information." << std::endl;
         return 1;
      }

   }
   catch (const std::runtime_error& ex)
   {
      Log("RunGfrd").fatal() << ex.what();
      return 2;
   }

   try
   {
      gfrd_print_header();
      Model m;
      switch (localDemoInit)
      {
         case 1:           // model file
         SetupModel(modelfile, model_arguments, world, rng, rules);
         break;

         case 2:
         {
            std::cout << std::setw(14) << "state file = " << modelfile << "\n";
            Persistence p;
            if (p.retreive(modelfile)) p.retreive_egfrd(*sptr);
            else THROW_EXCEPTION(std::runtime_error, "Could not loaded file.");
         } break;

         case 3:
         {
            std::cout << std::setw(14) << "crash file = " << modelfile << "\n";
            extSim.readSimFile(modelfile.c_str());     // try ASCII dump file (crash dump)
            if (extSim.active()) sptr = nullptr;            // render logic uses extSim when sptr == null!
            else THROW_EXCEPTION(std::runtime_error, "Could not loaded file.");
         } break;


         case 0:
         case 101: // Init GFRD (with three particle types in a box, cycling A -> B -> C -> A)
         {
            auto s1 = m.add_species_type(SpeciesType("A", m.get_def_structure_type_id(), 1e-12, 1e-9));
            auto s2 = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 1e-12, 0.5e-9));
            auto s3 = m.add_species_type(SpeciesType("C", m.get_def_structure_type_id(), 1e-12, 0.25e-9));

            world.initialize(1e-7, m);
            world.throwInParticles(s1, 24, rng, false);
            world.throwInParticles(s2, 14, rng, false);

            rules.add_reaction_rule(ReactionRule(s1, 1.0E-6, std::vector < SpeciesTypeID > {s2}));
            rules.add_reaction_rule(ReactionRule(s2, 1.0E-6, std::vector < SpeciesTypeID > {s3}));
            rules.add_reaction_rule(ReactionRule(s3, 1.0E-6, std::vector < SpeciesTypeID > {s1}));
         }
         break;

         case 102:            // demo DNA string, just for the fun
         {
            auto sG = m.add_structure_type(StructureType("Guanine"));
            auto sC = m.add_structure_type(StructureType("Cytosine"));
            auto sA = m.add_structure_type(StructureType("Adenine"));
            auto sT = m.add_structure_type(StructureType("Thymine"));
            auto sD = m.add_structure_type(StructureType("Helix"));
            std::map<StructureTypeID, StructureTypeID> complement = { {sG, sT}, {sT, sG}, {sA, sC} ,{sC, sA} };
            std::vector<StructureTypeID> sequence = { sG, sA, sT, sA, sA, sA, sT, sC, sT, sG, sG, sT, sC, sT, sT, sA, };

            auto sX = m.add_species_type(SpeciesType("X", m.get_def_structure_type_id(), 1e-12, 3e-9));

            world.initialize(1E-6, m);
            auto ws = world.world_size();
            auto wsid = world.get_def_structure_id();

            double turns = 0.75;
            double length = ws.X() / 8, radius = ws.X() / 80;
            for (size_t i = 0; i < sequence.size(); i++)
            {
               double f = static_cast<double>(i) / (sequence.size() - 1);
               Vector3 posY = Vector3(ws.X() / 2, f * ws.Y(), ws.Z() / 2);
               double phi = f * turns * 2 * M_PI;
               double arclength = std::sqrt(std::pow(M_PI * 2 * length*turns, 2) + std::pow(ws.Y(), 2)) / (2 * (sequence.size() - 1));
               Vector3 strand = length * Vector3(std::cos(phi), 0, std::sin(phi));
               Vector3 vz = Vector3::transformVector(Vector3::ux, Matrix4::createRotationY(phi));
               auto sID = sequence[i];
               double alpha = std::atan2(ws.Y() / turns, 2 * M_PI*length);

               world.add_structure(std::make_shared<CylindricalSurface>(CylindricalSurface("", sID, wsid, Cylinder(posY - 0.5*strand, radius, vz, length / 2))));
               world.add_structure(std::make_shared<CylindricalSurface>(CylindricalSurface("", complement[sID], wsid, Cylinder(posY + 0.5*strand, radius, vz, length / 2))));

               Vector3 vh = Vector3::transformVector(Vector3::uz, Matrix4::createRotationX(-alpha));
               Vector3 vhh1 = Vector3::transformVector(vh, Matrix4::createRotationY(phi));
               world.add_structure(std::make_shared<CylindricalSurface>(CylindricalSurface("Helix1", sD, wsid, Cylinder(posY - strand, 1.5 * radius, vhh1, arclength))));
               Vector3 vhh2 = Vector3::transformVector(vh, Matrix4::createRotationY(phi + M_PI));
               world.add_structure(std::make_shared<CylindricalSurface>(CylindricalSurface("Helix2", sD, wsid, Cylinder(posY + strand, 1.5 * radius, vhh2, arclength))));
            }

            world.throwInParticles(sX, 20, rng, false);
            //rules.add_interaction_rule(InteractionRule(sX, sD, 0.2, std::vector < SpeciesTypeID > {sX}));
         }
         break;

         case 103:        // Init GFRD (with PlanarSurfaces, not functional yet)
         {
            auto sPlane = m.add_structure_type(StructureType("plane"));

            auto sA = m.add_species_type(SpeciesType("A", sPlane, 2e-12, 1e-9));
            auto sB = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 1e-12, 3e-9));

            world.initialize(6e-6, m);
            auto ws = world.world_size();
            auto wsid = world.get_def_structure_id();

            auto vy = Vector3::transformVector(Vector3::uz, Matrix4::createRotationX(M_PI / 3.0));
            auto plane = std::make_shared<PlanarSurface>(PlanarSurface("plane", sPlane, wsid, Plane(ws / 2, Vector3::ux, vy, 0.3 * ws.X(), 0.3 * ws.Y(), false)));
            auto psid = world.add_structure(plane);
            UNUSED(psid);

            for (int i = 0; i < 3; ++i)
            {
               auto plane2 = std::make_shared<PlanarSurface>(PlanarSurface(make_string() << "plane" << i, sPlane, wsid, Plane(0.25 * ws + Vector3(0.2*i*ws.X(), 0, 0), Vector3::uz, Vector3::uy, 0.2 * ws.X(), 0.2*ws.Y(), false)));
               world.add_structure(plane2);
            }

            world.throwInParticles(sA, 2, rng, false, 0.2 * ws, 0.8 * ws);
            world.throwInParticles(sB, 1, rng, false);

            // Reaction Rules bind and unbind to plane
            rules.add_interaction_rule(InteractionRule(sB, sPlane, 0.2, std::vector < SpeciesTypeID > {sA}));
            rules.add_interaction_rule(InteractionRule(sA, world.get_def_structure_type_id(), 0.2, std::vector < SpeciesTypeID > {sB}));
         }
         break;

         case 104:
         {
            auto s1 = m.add_species_type(SpeciesType("A", m.get_def_structure_type_id(), 1e-12, 1e-9));
            auto s2 = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 1e-12, 0.5e-9));

            world.initialize(1e-7, m);

            // Test Multi construction
            world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(-3e-9, 0, 1e-8));
            world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(0, 1.5e-9, 1e-8));
            world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(3e-9, 0, 1e-8));
            world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(-6e-9, -1.5e-9, 1e-8));
            world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(6e-9, -1.5e-9, 1e-8));

            // Test Pair construction1
            const int particles = 8;
            for (int i = 0; i < particles; ++i)
            {
               world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3((-particles / 2 + i) * world.world_size().X() / (particles + 2), 0, -1e-8));
               world.add_particle(s2, world.get_def_structure_id(), world.world_size() / 2 + Vector3((-particles / 2 + i) * world.world_size().X() / (particles + 2), 4e-9, -1e-8));
            }
         } break;

          case 105:
          {
              // Test interaction between membrane-bound and 3D-diffusing particles


              auto sPlane = m.add_structure_type(StructureType("plane"));

              auto sA = m.add_species_type(SpeciesType("A", sPlane, 1e-12, 3e-8));
              auto sB = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 1e-11, 3e-8));
              auto sC = m.add_species_type(SpeciesType("C", sPlane, 5e-13, 3e-8));

              world.initialize(6e-6, m);
              auto ws = world.world_size();
              auto wsid = world.get_def_structure_id();

              auto origin = Vector3(0, 0, 0);
              auto pos = Vector3(ws.X() / 2, 0, ws.Z() / 2);
              auto uz = Vector3::uz;
//              auto plane = std::make_shared<PlanarSurface>(PlanarSurface("plane", sPlane, wsid, Plane(pos, Vector3::ux, uz, 0.5 * ws.X(), 0.5 * ws.Z(), false)));
//              auto psid = world.add_structure(plane);
//              UNUSED(psid);

//              world.throwInParticles(sA, 10, rng, false, Vector3(0, pos.Y(), 0), Vector3(ws.X(), pos.Y(), ws.Z()));
              world.throwInParticles(sB, 30, rng, false, Vector3(0, pos.Y(), 0), Vector3(ws.X(), ws.Y(), ws.Z()));

//              world.add_particle(sB, world.get_def_structure_id(), Vector3(4.2198775752446821e-06, 5.1948420752381338e-07, 5.9361989135308565e-06));
//              world.add_particle(sB, world.get_def_structure_id(), Vector3(4.1343342814910925e-06, 2.9702066521452007e-07, 5.9544249965700063e-06));

//            rules.add_reaction_rule(ReactionRule(sB, 0.21, std::vector<SpeciesTypeID>{sC}));
              rules.add_reaction_rule(ReactionRule(sA, sB, 0.21, std::vector<SpeciesTypeID>{sC}));
              rules.add_reaction_rule(ReactionRule(sB, sB, 0.21, std::vector<SpeciesTypeID>{sC}));

              // Reaction Rules bind and unbind to plane
//            rules.add_interaction_rule(InteractionRule(sB, sPlane, 0.2, std::vector < SpeciesTypeID > {sA}));
              rules.add_interaction_rule(InteractionRule(sA, world.get_def_structure_type_id(), 0.22, std::vector < SpeciesTypeID > {sB}));
          } break;

          case 106:
          {
              // Test PairSpherical domain bug occurring in Demo 5


              auto sA = m.add_species_type(SpeciesType("A", m.get_def_structure_type_id(), 1e-11, 3e-8));
              auto sB = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 1e-11, 3e-8));
              auto sC = m.add_species_type(SpeciesType("C", m.get_def_structure_type_id(), 1e-12, 1.5e-9));
              auto ws = world.world_size();

              world.initialize(6e-6, m);

              world.add_particle(sA, world.get_def_structure_id(), Vector3(4.2198775752446821e-06, 5.1948420752381338e-07, 5.9361989135308565e-06));
              world.add_particle(sB, world.get_def_structure_id(), Vector3(4.1343342814910925e-06, 2.9702066521452007e-07, 5.9544249965700063e-06));
              world.add_particle(sB, world.get_def_structure_id(), Vector3(4.569515338729998e-06, 5.7614426642224702e-07, 5.5921173272375657e-06));

              rules.add_reaction_rule(ReactionRule(sA, sB, 0.21, std::vector<SpeciesTypeID>{sC}));
          }
          break;

         default: THROW_EXCEPTION(illegal_size, "Demo number out of range.");

      }
   }
   catch (const std::runtime_error& ex)
   {
      Log("RunGfrd").fatal() << ex.what();
      return 2;
   }


   screen_capture_startup();

   // setup OpenGL
   glutInitWindowSize(800, 600);
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
   glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

   glutCreateWindow("GFRD Visualizer");

#if defined(_MSC_VER)
   HWND hWnd = FindWindowA(nullptr, "GFRD Visualizer");
   HANDLE hIcon = LoadIcon(GetModuleHandle(nullptr), MAKEINTRESOURCE(IDI_ICON1));
   SendMessage(hWnd, WM_SETICON, ICON_BIG, reinterpret_cast<LPARAM>(hIcon));
#endif

   glutMouseFunc(handleMouse);
   glutMotionFunc(handleMouseMotion);
   glutPassiveMotionFunc(handlePassiveMouseMotion);
   glutMouseWheelFunc(handleMouseWheel);
   glutKeyboardFunc(handleKeyboard);
   glutDisplayFunc(handleDisplay);
   glutReshapeFunc(handleReshape);
   glutIdleFunc(handleIdle);
   initRenderState();

   // rendering of -any- world size is scaled to GLunitlength (1.0) coordinate range
   auto ws = sptr ? sptr->world().world_size() : extSim.world_size();
   cam.set_distance(3);
   double side = 0.5 / std::max(ws.X(), std::max(ws.Y(), ws.Z()));      // scale largest side to GLunits 1.0
   cam.lookAt(side * ws);

   // run message loop
   glutMainLoop();

   screen_capture_shutdown();
   return 0;
}

// --------------------------------------------------------------------------------------------------------------------------------
