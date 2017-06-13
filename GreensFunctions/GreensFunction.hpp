#ifndef GREENSFUNCTION_HPP
#define GREENSFUNCTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "DefsGf.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT GreensFunction
{
public:
   enum class EventKind { IV_ESCAPE, IV_REACTION };

   GreensFunction(double D) noexcept : D_(D) {}
   virtual ~GreensFunction() = default;

   // prevent default construction, assignment and move operations, allow copy construction
   GreensFunction() = delete;
   GreensFunction(const GreensFunction&) = default;
   GreensFunction& operator=(const GreensFunction&) = default;
   GreensFunction(const GreensFunction&&) = delete;
   GreensFunction& operator=(const GreensFunction&&) = delete;

   double getD() const { return D_; }

   virtual std::string dump() const = 0;;
   virtual const char* type_name() const = 0;

   // common draw methods

   virtual double drawTime(double rnd) const = 0;

   virtual double drawR(double rnd, double t) const = 0;

   virtual EventKind drawEventType(double rnd, double t) const { UNUSED(rnd); UNUSED(t); return EventKind::IV_ESCAPE; }

protected:
   const double D_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif // GREENSFUNCTION_HPP