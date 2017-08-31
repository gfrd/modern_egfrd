#include <memory>

int main()
{
  auto check = std::make_unique<int>(42);
  return ((*check) == 42) ? 0 : 1;
}
