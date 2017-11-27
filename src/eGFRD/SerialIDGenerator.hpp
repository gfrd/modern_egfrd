#pragma once

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
