#ifndef MATRIX_SPACE_HPP
#define MATRIX_SPACE_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <array>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "exceptions.hpp"
#include "Vector3.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TObj, typename TKey, uint NX, uint NY, uint NZ>
class MatrixSpace
{
public:
   using key_type = TKey;
   using mapped_type = TObj;
   using value_type = std::pair<key_type, mapped_type>;
   using container_type = std::vector<value_type>;
   using size_type = uint;

private:
   using cell_type = std::vector<size_type>;
   using matrix_type = std::array<cell_type, NX*NY*NZ>;
   using cell_index_type = std::array<unsigned int, 3>;
   using cell_offset_type = std::array<signed int, 3>;
   using key_to_value_mapper_type = std::unordered_map<key_type, size_type>;

public:
   using iterator = typename container_type::iterator;
   using const_iterator = typename container_type::const_iterator;
   static const uint SizeX = NX;
   static const uint SizeY = NY;
   static const uint SizeZ = NZ;

   MatrixSpace() noexcept : world_size_(Vector3(1, 1, 1)), cell_size_(1.0 / NX) {}

   void initialize(double cell_size)
   {
      THROW_UNLESS_MSG(unsupported, size() == 0 || cell_size == cell_size_, "Cannot change size of non-empty matrix!.");
      cell_size_ = cell_size;
      world_size_ = Vector3(NX, NY, NZ) * cell_size;     // world is buildup with cubical cells (of cell_size) with NX, NY and NZ number of cells on each axis
   }

   // offset and matrix index calculation

   cell_index_type index(const Vector3& pos) const
   {
      return cell_index_type({ static_cast<unsigned int>(pos.X() / cell_size_) % NX,
         static_cast<unsigned int>(pos.Y() / cell_size_) % NY,
         static_cast<unsigned int>(pos.Z() / cell_size_) % NZ });
   }

private:

   bool offset_index(cell_index_type& i, const cell_offset_type& o) const
   {
      if ((o[0] < 0 && i[0] + o[0] < 0) || (i[0] + o[0] >= NX) ||
         (o[1] < 0 && i[1] + o[1] < 0) || (i[1] + o[1] >= NY) ||
         (o[2] < 0 && i[2] + o[2] < 0) || (i[2] + o[2] >= NZ))
      {
         return false;
      }
      i[0] += o[0];
      i[1] += o[1];
      i[2] += o[2];
      return true;
   }

   Vector3 offset_index_cyclic(cell_index_type& i, const cell_offset_type& o) const
   {
      double x = 0, y = 0, z = 0;
      if (o[0] < 0 && static_cast<int>(i[0] + o[0]) < 0)
      {
         size_type t = (i[0] + NX - (-o[0] % NX)) % NX;
         x = (o[0] - static_cast<int>(t - i[0])) * cell_size_;
         i[0] = t;
      }
      else if (i[0] + o[0] >= NX)
      {
         size_type t = (i[0] + (o[0] % NX)) % NX;
         x = (o[0] - (t - i[0])) * cell_size_;
         i[0] = t;
      }
      else
         i[0] += o[0];

      if (o[1] < 0 && static_cast<int>(i[1] + o[1]) < 0)
      {
         size_type t = (i[1] + NY - (-o[1] % NY)) % NY;
         y = (o[1] - static_cast<int>(t - i[1])) * cell_size_;
         i[1] = t;
      }
      else if (i[1] + o[1] >= NY)
      {
         size_type t = (i[1] + (o[1] % NY)) % NY;
         y = (o[1] - (t - i[1])) * cell_size_;
         i[1] = t;
      }
      else
         i[1] += o[1];

      if (o[2] < 0 && static_cast<int>(i[2] + o[2]) < 0)
      {
         size_type t = (i[2] + NZ - (-o[2] % NZ)) % NZ;
         z = (o[2] - static_cast<int>(t - i[2])) * cell_size_;
         i[2] = t;
      }
      else if (i[2] + o[2] >= NZ)
      {
         size_type t = (i[2] + (o[2] % NZ)) % NZ;
         z = (o[2] - (t - i[2])) * cell_size_;
         i[2] = t;
      }
      else
         i[2] += o[2];

      return Vector3(x, y, z);
   }

public:

   // Get basic properties

   const cell_type& cell(const cell_index_type& i) const { return matrix_[i[0] + i[1] * NX + i[2] * NX*NY]; }

   cell_type& cell(const cell_index_type& i) { return matrix_[i[0] + i[1] * NX + i[2] * NX*NY]; }

   Vector3 world_size() const { return world_size_; }

   double cell_size() const { return cell_size_; }

   std::array<size_type, 3> matrix_size() const { return std::array<size_type, 3>({ NX,NY,NZ }); }

   size_type size() const { return static_cast<size_type>(values_.size()); }

   iterator update(const iterator& old_value, const value_type& v)
   {
      const Vector3& pos = v.second.position();
      ASSERT(0.0 <= pos.X() && pos.X() < world_size_.X());
      ASSERT(0.0 <= pos.Y() && pos.Y() < world_size_.Y());   // When one of these fire, you try to insert a position that does not fit inside the matrix!
      ASSERT(0.0 <= pos.Z() && pos.Z() < world_size_.Z());

      cell_type* new_cell = &cell(index(pos));
      cell_type* old_cell = nullptr;

      if (old_value != values_.end()) old_cell = &cell(index((*old_value).second.position()));

      if (new_cell == old_cell)
      {
         *old_value = v;
         return old_value;
      }

      size_type index;
      if (old_cell)
      {
         *old_value = v;

         auto i = std::find(old_cell->begin(), old_cell->end(), (old_value - values_.begin()));
         index = *i;
         old_cell->erase(i);
      }
      else
      {
         index = size();
         values_.emplace_back(v);
         rmap_[v.first] = index;
      }
      new_cell->emplace_back(index); //new_cell->push(index);
      return values_.begin() + index;
   }

   std::pair<iterator, bool> update(const value_type& v)
   {
      cell_type* new_cell = &cell(index(v.second.position()));
      auto &&old_value = values_.end();
      cell_type* old_cell = nullptr;

      {
         auto &&i = rmap_.find(v.first);
         if (i != rmap_.end())
         {
            old_value = values_.begin() + (*i).second;
            old_cell = &cell(index(old_value->second.position()));
         }
      }

      if (new_cell == old_cell)
      {
         *old_value = v;
         return std::pair<iterator, bool>(old_value, false);
      }

      size_type index;
      if (old_cell)
      {
         *old_value = v;

         auto i = std::find(old_cell->begin(), old_cell->end(), (old_value - values_.begin()));
         index = *i;
         old_cell->erase(i);
      }
      else
      {
         index = size();
         values_.emplace_back(v);
         rmap_[v.first] = index;
      }

      new_cell->emplace_back(index);
      return std::pair<iterator, bool>(values_.begin() + index, !old_cell);
   }

   bool erase(const iterator& i)
   {
      if (end() == i) return false;
      const size_type old_index = static_cast<size_type>(i - values_.begin());

      cell_type& e_cell = cell(index((*i).second.position()));
      e_cell.erase(std::find(e_cell.begin(), e_cell.end(), old_index));       //e_cell.erase(old_index);
      rmap_.erase((*i).first);

      size_type const last_index = size() - 1;
      if (old_index < last_index)
      {
         const value_type& last = values_[last_index];
         cell_type& old_c = cell(index(last.second.position()));
         old_c.erase(std::find(old_c.begin(), old_c.end(), last_index));       //old_c.erase(last_index);
         old_c.emplace_back(old_index);  //old_c.push(old_index);
         rmap_[last.first] = old_index;
         *i = last;
      }
      values_.pop_back();
      return true;
   }

   bool erase(const key_type& k)
   {
      auto &&p = rmap_.find(k);
      if (rmap_.end() == p) return false;
      return erase(values_.begin() + (*p).second);
   }

   void clear()
   {
      for (auto &p : matrix_)
         p.clear();
      rmap_.clear();
      values_.clear();
   }

   // Iterations

   iterator begin() { return values_.begin(); }

   const_iterator begin() const { return values_.cbegin(); }

   iterator end() { return values_.end(); }

   const_iterator end() const { return values_.cend(); }

   // Find

   iterator find(const key_type& k)
   {
      auto &&p = rmap_.find(k);
      if (rmap_.end() == p)
         return values_.end();
      return values_.begin() + (*p).second;
   }

   const_iterator find(const key_type& k) const
   {
      auto &&p = rmap_.find(k);
      if (rmap_.end() == p)
         return values_.end();
      return values_.begin() + (*p).second;
   }

   // Neighbor(s)

   template<typename TCollect>
   void each_neighbor(const cell_index_type& idx, TCollect& collector) const
   {
      cell_offset_type off;
      for (off[2] = -1; off[2] <= 1; ++off[2])
      {
         for (off[1] = -1; off[1] <= 1; ++off[1])
         {
            for (off[0] = -1; off[0] <= 1; ++off[0])
            {
               cell_index_type _idx = idx;
               if (!offset_index(_idx, off)) continue;
               const cell_type& c = cell(_idx);
               for (const auto &i : c) collector(values_.cbegin() + i, Vector3());
            }
         }
      }
   }

   template<typename TCollect>
   void each_neighbor_cyclic(const cell_index_type& idx, TCollect& collector) const
   {
      cell_offset_type off;
      for (off[2] = -1; off[2] <= 1; ++off[2])
      {
         for (off[1] = -1; off[1] <= 1; ++off[1])
         {
            for (off[0] = -1; off[0] <= 1; ++off[0])
            {
               cell_index_type _idx = idx;
               const Vector3 pos_off = offset_index_cyclic(_idx, off);
               const cell_type& c = cell(_idx);
               for (const auto &i : c) collector(values_.cbegin() + i, pos_off);
            }
         }
      }
   }

private:
   friend class Persistence;

   Vector3 world_size_;
   double cell_size_;
   matrix_type matrix_;
   key_to_value_mapper_type rmap_;
   container_type values_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* MATRIX_SPACE_HPP */
