#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <GL/freeglut.h>
#include "Vector2.hpp"
#include "Vector3.hpp"
#include "Matrix4.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

const double mouseDragSpeed = 6.5;

// --------------------------------------------------------------------------------------------------------------------------------

class CameraController
{
public:
   CameraController() noexcept : theta_(0), azimuth_(0), distance_(10), lookAt_(Vector3()), mouseDrag_(false),
      downX_(0), downY_(0), value1_(0), value2_(0) {}
   virtual ~CameraController() = default;

   // glut entry functions
   void handleMouseWheel(int wheel, int direction, int x, int y)
   {
      UNUSED(wheel, x, y);
      distance_ += static_cast<float>(0.05 * direction * (distance_ == 0 ? 0.1 : distance_));
      if (distance_ < 0.0) distance_ = 0.0;
      glutPostRedisplay();
   }

   virtual void handleMouse(int button, int updown, int x, int y)
   {
      UNUSED(button);
      if (updown == GLUT_DOWN && !mouseDrag_)
      {
         mouseDrag_ = true;
         downX_ = x;
         downY_ = y;
         value1_ = azimuth_;
         value2_ = theta_;
      }

      if (updown == GLUT_UP && mouseDrag_)
      {
         mouseDrag_ = false;
         glutPostRedisplay();
      }
   }

   virtual void handleMouseMotion(int x, int y)
   {
      if (mouseDrag_)
      {
         azimuth_ = value1_ + (x - downX_) / (100 * mouseDragSpeed);
         theta_ = value2_ + (y - downY_) / (100 * mouseDragSpeed);
         glutPostRedisplay();
      }
   }

   virtual void handlePasiveMouseMotion(int x, int y)
   {
      UNUSED(x, y);
      glutSetCursor(GLUT_CURSOR_CROSSHAIR);
   }

   double distance() const { return distance_; }
   double theta() const { return theta_; }
   double azimuth() const { return azimuth_; }

   void set_distance(double distance) { distance_ = std::abs(distance); }
   void set_angles(double azimuth, double theta) { azimuth_ = azimuth; theta_ = theta; }

   Vector3 eye() const { return position() + center(); }

   const Vector3& center() const { return lookAt_; }

   void virtual lookAt(const Vector3& look) { lookAt_ = look; }

   static const Vector3& Up() { return Vector3::uy; }

   bool drag() const { return mouseDrag_; }

   Vector3 position() const
   {
      return calculate(Vector3(0, 0, distance_));
   }

   Vector3 calculate(const Vector3& u) const
   {
      auto v1 = Vector3::transformVector(u, Matrix4::createRotationX(theta_));
      return Vector3::transformVector(v1, Matrix4::createRotationY(azimuth_));
   }

   void look() const
   {
      auto e = eye();
      gluLookAt(e.X(), e.Y(), e.Z(), lookAt_.X(), lookAt_.Y(), lookAt_.Z(), 0, 1, 0);
   }

   void virtual lookTo(const Vector3& lookTo)
   {
      auto p2 = eye() - lookTo;

      double distance = p2.length();
      double azimuth = std::atan2(-p2.X(), p2.Z());
      auto v1 = Vector3::transformVector(p2, Matrix4::createRotationY(-azimuth));
      double theta = std::atan2(v1.Y(), v1.Z());

      distance_ = distance;
      theta_ = theta;
      azimuth_ = azimuth;
      lookAt_ = lookTo;
   }

protected:
   double theta_;
   double azimuth_;
   double distance_;
   Vector3 lookAt_;
   bool mouseDrag_;
   int downX_, downY_;
   double value1_;
   double value2_;
};

// --------------------------------------------------------------------------------------------------------------------------------

struct Rect
{
   Rect() noexcept : x_(0), y_(0), width_(0), height_(0) {}
   Rect(int x, int y, int width, int height) noexcept : x_(x), y_(y), width_(width), height_(height) {}

   bool contains(int x, int y) const { return x >= x_ && x < x_ + width_ && y >= y_ && y < y_ + height_; }

private:
   int x_, y_, width_, height_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class CameraControllerExtended : public CameraController
{
public:
   CameraControllerExtended() noexcept : CameraController(), dragMode_(DragMode::None) { }

   // glut entry functions
   void handleReshape(int cx, int cy)
   {
      // Calculate Drag area's
      areaRotate_ = Rect(cx / 10, cy / 10, 8 * cx / 10, 8 * cy / 10);	      // Middle area
      areaDepth_ = Rect(cx / 10, 0, 8 * cx / 10, cy / 10);					   // Middle top side of window;
      areaShiftX_ = Rect(cx / 10, 9 * cy / 10, 8 * cx / 10, cy / 10);		   // Middle bottom side of window;
      areaShiftZ_ = Rect(0, 0, cx / 10, cy);										   // left side of window
      areaShiftY_ = Rect(9 * cx / 10, 0, cx / 10, cy);						      // right side of window
   }

   virtual void handleMouse(int button, int updown, int x, int y) override
   {
      CameraController::handleMouse(button, updown, x, y);
      if (mouseDrag_)
      {
         if (areaRotate_.contains(x, y))
         {
            dragMode_ = DragMode::Rotate;
            value1_ = azimuth_;
            value2_ = theta_;
         }
         else if (areaDepth_.contains(x, y))
         {
            dragMode_ = DragMode::Depth;
            value1_ = distance_;
         }
         else if (areaShiftX_.contains(x, y))
         {
            dragMode_ = DragMode::ShiftX;
            value1_ = lookAt_.X();
            value2_ = lookAt_.Z();
         }
         else if (areaShiftY_.contains(x, y))
         {
            dragMode_ = DragMode::ShiftY;
            value1_ = lookAt_.Y();
         }
         else if (areaShiftZ_.contains(x, y))
         {
            dragMode_ = DragMode::ShiftZ;
            value1_ = lookAt_.X();
            value2_ = lookAt_.Z();
         }
      }
      else dragMode_ = DragMode::None;
   }

   virtual void handleMouseMotion(int x, int y) override
   {
      switch (dragMode_)
      {
      default:
      case DragMode::None:
         break;

      case DragMode::Rotate:
         azimuth_ = value1_ + (x - downX_) / (100 * mouseDragSpeed);
         theta_ = value2_ + (y - downY_) / (100 * mouseDragSpeed);
         glutPostRedisplay();
         break;

      case DragMode::Depth:
         distance_ = value1_ + (downX_ - x) / (10 * mouseDragSpeed);
         if (distance_ < 0.0) distance_ = 0.0;
         glutPostRedisplay();
         break;

      case DragMode::ShiftX:
         lookAt(Vector3(value1_ + (downX_ - x) / (60 * mouseDragSpeed) * std::cos(azimuth_), lookAt_.Y(), value2_ + (downX_ - x) / (60 * mouseDragSpeed) * std::sin(azimuth_)));
         glutPostRedisplay();
         break;

      case DragMode::ShiftY:
         lookAt(Vector3(lookAt_.X(), value1_ - ((downY_ - y) / (60 * mouseDragSpeed)), lookAt_.Z()));
         glutPostRedisplay();

         break;

      case DragMode::ShiftZ:
         lookAt(Vector3(value1_ + (downY_ - y) / (60 * mouseDragSpeed) * std::sin(azimuth_), lookAt_.Y(), value2_ - (downY_ - y) / (60 * mouseDragSpeed) * std::cos(azimuth_)));
         glutPostRedisplay();
         break;
      }
   }

   virtual void handlePasiveMouseMotion(int x, int y) override
   {
      if (areaRotate_.contains(x, y))
         glutSetCursor(GLUT_CURSOR_CROSSHAIR);
      else if (areaDepth_.contains(x, y))
         glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
      else if (areaShiftX_.contains(x, y))
         glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
      else if (areaShiftY_.contains(x, y))
         glutSetCursor(GLUT_CURSOR_UP_DOWN);
      else if (areaShiftZ_.contains(x, y))
         glutSetCursor(GLUT_CURSOR_UP_DOWN);
   }

protected:
   enum class DragMode { None = 0, Rotate, Depth, ShiftX, ShiftY, ShiftZ };

   DragMode dragMode_;
   Rect areaRotate_;
   Rect areaDepth_;
   Rect areaShiftX_;
   Rect areaShiftY_;
   Rect areaShiftZ_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class CameraControllerAnimated : public CameraControllerExtended
{
   const double DAMPING = 0.75;

public:
   CameraControllerAnimated() : CameraControllerExtended(), last1_(0), last2_(0), demo_(false), lookTo_(lookAt_) {}

   virtual void handleMouse(int button, int updown, int x, int y) override
   {
      if (!mouseDrag_)
      {
         speed_ = Vector2();
         last2_ = clock();
         pos_ = Vector2(azimuth_, theta_);
      }
      if (mouseDrag_ && dragMode_ == DragMode::Rotate)
      {
         CalcSpeed();
         anim_ = speed_;
         last1_ = clock();
      }
      CameraControllerExtended::handleMouse(button, updown, x, y);
   }

   virtual void handleMouseMotion(int x, int y) override
   {
      CameraControllerExtended::handleMouseMotion(x, y);
      if (mouseDrag_ && dragMode_ == DragMode::Rotate)
         CalcSpeed();
   }

   virtual void handleIdle()
   {
      if (mouseDrag_ || ((anim_.length() < 1E-3) && ((lookTo_ - lookAt_).length() < 1E-3))) return;

      clock_t time = clock();
      double dt = static_cast<double>(time - last1_) / CLOCKS_PER_SEC;
      last1_ = time;

      if ((lookTo_ - lookAt_).length() > 1E-3) lookTo(lookTo_, dt);
      //else
      if (anim_.length() > 1E-3)
      {
         if (demo_) anim_ = Vector2(0.5, 0.25 * std::sin(time / 2000));

         anim_ -= anim_ * (DAMPING * dt);
         azimuth_ += dt * anim_.X();
         theta_ += dt * anim_.Y();
      }

      glutPostRedisplay();
   }

   void set_demo(bool demo) { demo_ = demo; if (demo_) { anim_ = Vector2(1.0, 0.0); last1_ = clock(); } }
   bool demo() const { return demo_; }

   void lookAt(const Vector3& look) override { lookAt_ = look; lookTo_ = look; }

   void lookTo(const Vector3& lookTo) override { lookTo_ = lookTo; last1_ = clock(); }

private:
   void CalcSpeed()
   {
      clock_t time = clock();
      if (time == last2_) return;   // skip zero time
      double dt = static_cast<double>(time - last2_) / CLOCKS_PER_SEC;
      last2_ = time;

      Vector2 pos(Vector2(azimuth_, theta_));
      speed_ = (pos - pos_) / dt;
      pos_ = pos;
   }

   void lookTo(const Vector3& lookTo, double dt)
   {
      auto e = eye();
      auto p2 = e - lookTo;
      double distance = p2.length();
      double azimuth = std::atan2(-p2.X(), p2.Z());
      auto v1 = Vector3::transformVector(p2, Matrix4::createRotationY(-azimuth));
      double theta = std::atan2(v1.Y(), v1.Z());

      distance_ += (distance - distance_) * (DAMPING * dt);
      theta_ += shortAngle(theta , theta_) * (DAMPING * dt);
      azimuth_ += shortAngle(azimuth , azimuth_) * (DAMPING * dt);
      lookAt_ = e - position();
   }

   double shortAngle(double x, double y) const
   {
      return std::atan2(std::sin(x - y), std::cos(x - y));
   }

protected:
   Vector2 anim_;
   Vector2 speed_;
   Vector2 pos_;
   clock_t last1_;
   clock_t last2_;
   bool demo_;
   Vector3 lookTo_;
};

// --------------------------------------------------------------------------------------------------------------------------------
