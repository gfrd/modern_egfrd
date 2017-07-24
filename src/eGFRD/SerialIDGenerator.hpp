#ifndef SERIAL_ID_GENERATOR_HPP
#define SERIAL_ID_GENERATOR_HPP

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TIdentifier>
struct SerialIDGenerator
{
   SerialIDGenerator() : next_(TIdentifier(1u)) {}
   TIdentifier operator()() { return TIdentifier(++next_); }

   friend class Persistence;

private:
   TIdentifier next_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* SERIAL_ID_GENERATOR_HPP */
