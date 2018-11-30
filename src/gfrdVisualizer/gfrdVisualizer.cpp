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
std::string statefile;

// --------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> help = { "ESC = exit" , "H = show help", "F = full screen", "S = show shells", "D = demo rotate", "P = screenshot", "O = draw origin" ,"I = particle ID/SID" ,
                                   ". = sim step", "/ = check sim", "A = auto", "B = burst all",  "</> = select structure",
                                   "+/- = select particle", "g = goto particle", "0 = deselect all", "[/] = select domain", "j = goto domain", "C = center world", };

// --------------------------------------------------------------------------------------------------------------------------------

size_t get_help_width()
{
   size_t width = 0;
   for (auto& txt : help)
   {
      size_t size = glutBitmapWidth(GLUT_BITMAP_9_BY_15, 'c') * txt.length();
      width = std::max(width, size);
   };
   return width;
}

// --------------------------------------------------------------------------------------------------------------------------------

void check_sim()
{
   try
   {
      sptr->check();
   }
   catch (gfrd_exception ex)
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

      case '.': if (sptr) sptr->step(); else extSim.readSimFile(extSim.get_filename()); glutPostRedisplay(); break;
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

      //case 'z': case 'Z': { Persistence p; p.store("d:\\simstate.bin"); p.store_egfrd(*sptr); } break;
      case 'x': case 'X': if (!statefile.empty()) { Persistence p; p.retreive(statefile); p.retreive_egfrd(*sptr); glutPostRedisplay(); } break;

      default: break;
      }
   }
   catch (std::exception ex)
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

void handleIdle(void)
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
   catch (std::exception ex)
   {
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


void drawOrigin(void)
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

void handleDisplay(void)
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
               glRasterPos2i(viewport[2] - static_cast<GLint>(help_width), y - glutBitmapHeight(GLUT_BITMAP_9_BY_15) + viewport[3]);
               glutBitmapString(GLUT_BITMAP_9_BY_15, reinterpret_cast<const unsigned char*>(txt.c_str()));
               y -= 18;
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
   catch (std::exception ex)
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

void handlePasiveMouseMotion(int x, int y)
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
   screen_capture_startup();
   glutInit(&argc, argv);
   glutInitWindowSize(800, 600);

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
   glutPassiveMotionFunc(handlePasiveMouseMotion);
   glutMouseWheelFunc(handleMouseWheel);
   glutKeyboardFunc(handleKeyboard);
   glutDisplayFunc(handleDisplay);
   glutReshapeFunc(handleReshape);
   glutIdleFunc(handleIdle);

   initRenderState();

   // local to main, need to stay in scope
   RandomNumberGenerator rng;
   ReactionRuleCollection rules;
   World world;
   EGFRDSimulator s(world, rules, rng);
   sptr = &s;

   bool localInit = true;
   if (argc == 2)
   {
      // try ASCI visual dump
      extSim.readSimFile(argv[1]);

      // try binary simulator dump
      if (extSim.active())
      {
         localInit = false;
         sptr = nullptr;      // render logic uses extSim when sptr == null!
      }
      else
         try
      {
         Persistence p;
         if (p.retreive(argv[1]))
         {
            p.retreive_egfrd(s);
            localInit = false;
            statefile = argv[1];
         }
      }
      catch (std::runtime_error ex)
      {
         Log("Visualize").warn() << "Failed to restore simulator state from file: '" << argv[1] << "' fault: " << ex.what();
      }
   }

   if (localInit)
   {
      switch (2)      // select demo init
      {
      default:
      case 0:        // Init GFRD (with PlanarSurfaces, not functional yet)
      {
         Model m;
         auto sPlane = m.add_structure_type(StructureType("plane"));

         auto sA = m.add_species_type(SpeciesType("A", sPlane, 2e-12, 1e-9));
         auto sB = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 1e-12, 3e-9));

         world.initialize(1E-6, m);
         auto ws = world.world_size();
         auto wsid = world.get_def_structure_id();

         auto vy = Vector3::transformVector(Vector3::uz, Matrix4::createRotationX(M_PI / 3.0));
         auto plane = std::make_shared<PlanarSurface>(PlanarSurface("plane", sPlane, wsid, Plane(ws / 2, Vector3::ux, vy, 0.3 * ws.X(), 0.3 * ws.Y(), false)));
         auto psid = world.add_structure(plane);
         UNUSED(psid);

         //for (int i = 0; i < 3; ++i)
         //{
         //   auto plane2 = std::make_shared<PlanarSurface>(PlanarSurface(make_string() << "plane" << i, sPlane, wsid, Plane(0.25 * ws + Vector3(0.2*i*ws.X(), 0, 0), Vector3::uz, Vector3::uy, 0.2 * ws.X(), 0.2*ws.Y(), false)));
         //   world.add_structure(plane2);
         //}

         world.throwInParticles(sA, 2, rng, false, 0.2 * ws, 0.8 * ws);
         world.throwInParticles(sB, 1, rng, false);

         // Reaction Rules bind and unbind to plane
         rules.add_interaction_rule(InteractionRule(sB, sPlane, 0.2, std::vector < SpeciesTypeID > {sA}));
         rules.add_interaction_rule(InteractionRule(sA, world.get_def_structure_type_id(), 0.2, std::vector < SpeciesTypeID > {sB}));
      }
      break;

      case 1:            // demo DNA string, just for the fun
      {
         Model m;
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


      case 2: // Init GFRD (with three particle types in a box, cycling A -> B -> C -> A)
      {
         Model m;
         auto s1 = m.add_species_type(SpeciesType("A", m.get_def_structure_type_id(), 1e-12, 1e-9));
         auto s2 = m.add_species_type(SpeciesType("B", m.get_def_structure_type_id(), 1e-12, 0.5e-9));
         auto s3 = m.add_species_type(SpeciesType("C", m.get_def_structure_type_id(), 1e-12, 0.25e-9));

         world.initialize(1e-7, m);

         world.throwInParticles(s1, 24, rng, false);
         world.throwInParticles(s2, 14, rng, false);

         // Test Multi construction
         //world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(-3e-9, 0, 0));
         //world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(0, 1.5e-9, 0));
         //world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(3e-9, 0, 0));
         //world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(-6e-9, -1.5e-9, 0));
         //world.add_particle(s1, world.get_def_structure_id(), world.world_size() / 2 + Vector3(6e-9, -1.5e-9, 0));

         // Test Pair construction
         //const int particles = 9;
         //for (int i = 0; i < particles; ++i)
         //{
         //    w.add_particle(s1, w.get_def_structure_id(), w.world_size() / 2 + Vector3((-particles / 2 + i) * w.world_size().X() / (particles + 2), 0, 0));
         //    w.add_particle(s2, w.get_def_structure_id(), w.world_size() / 2 + Vector3((-particles / 2 + i) * w.world_size().X() / (particles + 2), 4e-9, 0));
         //}

         //// Test Pair construction2
         //double wx2 = w.world_size().X() / 2;
         //w.add_particle(s1, w.get_def_structure_id(), w.world_size() / 2 + Vector3(wx2 - 1.1e-9, 0, 0));
         //w.add_particle(s1, w.get_def_structure_id(), w.world_size() / 2 - Vector3(wx2 - 1.1e-9, 0, 0));

         rules.add_reaction_rule(ReactionRule(s1, 1.0E-6, std::vector < SpeciesTypeID > {s2}));
         rules.add_reaction_rule(ReactionRule(s2, 1.0E-6, std::vector < SpeciesTypeID > {s3}));
         rules.add_reaction_rule(ReactionRule(s3, 1.0E-6, std::vector < SpeciesTypeID > {s1}));
      }
      break;
      }
   }

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
