#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <utility>
#include <vector>
#include <unordered_map>
#include "DefsEgfrd.hpp"
#include "SerialIDGenerator.hpp"
#include "makeString.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
class DynamicPriorityQueue
{
public:

   bool empty() const { return items_.empty(); }

   size_t size() const { return items_.size(); }

   std::pair<TID, const TItem&> top() const { return items_[top_index()]; }

   std::pair<TID, const TItem&> second() const { return items_[second_index()]; }

   const TItem& get(TID id) const { return items_[findIndex(id)].second; }

   void pop() { pop_by_index(top_index()); }

   void pop(TID id) { pop_by_index(findIndex(id)); }

   void clear()
   {
      items_.clear();
      heap_.clear();
      position_vector_.clear();
      index_map_.clear();
   }

   void replace(const TID id, TItem& value);

   TID push(const TItem& item);

   bool check() const;

   bool has(const TID id) const { return index_map_.find(id) != index_map_.end(); }

   //-----------------------------------------------------------------------------------------------------------

protected:
   size_t top_index() const
   {
      THROW_UNLESS(std::out_of_range, size() > 0);
      return heap_[0];
   }

   size_t second_index() const
   {
      THROW_UNLESS(std::out_of_range, size() > 1);

      size_t index1(heap_[1]);
      if (size() == 2) return index1;

      size_t index2(heap_[2]);
      return compare_(items_[index1].second, items_[index2].second) ? index1 : index2;
   }

   void pop_by_index(size_t index);

   void move_pos(size_t pos);

   void move_up_pos(size_t position, size_t start = 0);
   void move_up_pos_impl(size_t position, size_t start = 0);

   void move_down_pos(size_t position);
   void move_down_pos_impl(size_t position);

   bool check_size() const;
   bool check_position_mapping() const;
   bool check_heap() const;

   //-----------------------------------------------------------------------------------------------------------

   size_t findIndex(const TID& id) const
   {
      auto i(index_map_.find(id));
      if (i == index_map_.end()) throw std::out_of_range(make_string() << "Item with " << id << " not found.");
      return (*i).second;
   }

   TID storeIndex(size_t index)
   {
      const auto id(idgen_());
      index_map_.insert(typename std::unordered_map<TID, size_t>::value_type(id, index));
      return id;
   }

   void freeIndex(size_t index, TID id, TID last_item_id)
   {
      index_map_[last_item_id] = index;
      index_map_.erase(id);
   }

private:
   friend class Persistence;

   std::vector<std::pair<TID, TItem>> items_;
   std::vector<size_t> heap_;
   std::vector<size_t> position_vector_;
   std::unordered_map<TID, size_t> index_map_;
   SerialIDGenerator<TID> idgen_;
   typename TItem::comparator compare_;
};

// --------------------------------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
void DynamicPriorityQueue<TID, TItem>::move_pos(size_t pos)
{
   size_t index(heap_[pos]);
   const std::pair<TID, TItem>& item(items_[index]);
   size_t succ(2 * pos + 1);
   if (succ < size())
   {
      if (compare_(items_[heap_[succ]].second, item.second) || (succ + 1 < size() && compare_(items_[heap_[succ + 1]].second, item.second)))
      {
         move_down_pos_impl(pos);
         return;
      }
   }

   move_up_pos(pos);
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
void DynamicPriorityQueue<TID, TItem>::move_up_pos(size_t position, size_t start)
{
   if (position == 0) return;

   size_t index(heap_[position]);
   const std::pair<TID, TItem>& item(items_[index]);

   size_t pred((position - 1) / 2);
   size_t predindex_type(heap_[pred]);

   if (compare_(item.second, items_[predindex_type].second))
      move_up_pos_impl(position, start);
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
void DynamicPriorityQueue<TID, TItem>::move_down_pos(size_t position)
{
   size_t index(heap_[position]);
   const std::pair<TID, TItem>& item(items_[index]);

   size_t succ(2 * position + 1);
   if (succ < size())
      if (compare_(items_[heap_[succ]].second, item.second) || (succ + 1 < size() && compare_(items_[heap_[succ + 1]].second, item.second)))
         move_down_pos_impl(position);
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
void DynamicPriorityQueue<TID, TItem>::move_up_pos_impl(size_t position, size_t start)
{
   size_t index(heap_[position]);
   const std::pair<TID, TItem>& item(items_[index]);

   if (position <= start) return;

   size_t pos(position);
   size_t pred((pos - 1) / 2);
   size_t predindex_type(heap_[pred]);

   do
   {
      heap_[pos] = predindex_type;
      position_vector_[predindex_type] = pos;
      pos = pred;

      if (pos <= start) break;

      pred = (pos - 1) / 2;
      predindex_type = heap_[pred];

   } while (!compare_(items_[predindex_type].second, item.second));

   heap_[pos] = index;
   position_vector_[index] = pos;
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
void DynamicPriorityQueue<TID, TItem>::move_down_pos_impl(size_t position)
{
   size_t index(heap_[position]);
   size_t succ(2 * position + 1);
   size_t pos(position);

   while (succ < size())
   {
      size_t right_pos(succ + 1);
      if (right_pos < size() && !compare_(items_[heap_[succ]].second, items_[heap_[right_pos]].second))
         succ = right_pos;

      heap_[pos] = heap_[succ];
      position_vector_[heap_[pos]] = pos;
      pos = succ;
      succ = 2 * pos + 1;
   }

   heap_[pos] = index;
   position_vector_[index] = pos;

   move_up_pos(pos, position);
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
TID DynamicPriorityQueue<TID, TItem>::push(const TItem& item)
{
   size_t index(items_.size());
   TID id(storeIndex(index));
   items_.emplace_back(std::pair<TID, TItem>(id, item));
   heap_.emplace_back(index);
   position_vector_.emplace_back(index);
   move_up_pos(index);
   return id;
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
void DynamicPriorityQueue<TID, TItem>::pop_by_index(size_t index)
{
   std::pair<TID, TItem>& item(items_[index]);
   // 1. update index<->identifier_type mapping.
   freeIndex(index, item.first, items_.back().first);

   // 2. pop the item from the items_.
   std::swap(item, items_.back());
   items_.pop_back();

   size_t removed_pos(position_vector_[index]);
   size_t moved_pos(position_vector_.back());

   // 3. swap position_vector_[end] and position_vector_[index]
   position_vector_[index] = moved_pos;
   heap_[moved_pos] = index;

   // 4. if heap_[end] and heap_[removed] do not overlap, swap these, pop back, and update the heap_.
   if (removed_pos != heap_.size() - 1)
   {
      heap_[removed_pos] = heap_.back();
      position_vector_[heap_.back()] = removed_pos;

      position_vector_.pop_back();
      heap_.pop_back();

      move_pos(removed_pos);
   }
   else  // if heap_[end] and heap_[removed] are the same, simply pop back.
   {
      position_vector_.pop_back();
      heap_.pop_back();
   }
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
void DynamicPriorityQueue<TID, TItem>::replace(const TID id, TItem& value)
{
   size_t index(findIndex(id));
   items_[index].second = value;
   move_pos(position_vector_[index]);
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
bool DynamicPriorityQueue<TID, TItem>::check() const
{
   bool result(true);
   result = result && check_size();
   result = result && check_position_mapping();
   result = result && check_heap();
   return result;
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
bool DynamicPriorityQueue<TID, TItem>::check_size() const
{
   bool result(true);
   result = result && items_.size() == size();
   result = result && heap_.size() == size();
   result = result && position_vector_.size() == size();
   return result;
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
bool DynamicPriorityQueue<TID, TItem>::check_position_mapping() const
{
   bool result(true);
   for (size_t i(0); i < size(); ++i)
   {
      result = result && heap_[i] < size();
      result = result && position_vector_[i] < size();
      result = result && heap_[position_vector_[i]] == i;
   }
   return result;
}

//-----------------------------------------------------------------------------------------------------------

template<typename TID, typename TItem>
bool DynamicPriorityQueue<TID, TItem>::check_heap() const
{
   bool result(true);

   // assert correct ordering of items in the heap_.
   for (size_t pos(0); pos < size(); ++pos)
   {
      const std::pair<TID, TItem>& item(items_[heap_[pos]]);
      size_t succ(pos * 2 + 1);
      if (succ < size())
      {
         result = result && compare_(item.second, items_[heap_[succ]].second);
         size_t right_pos(succ + 1);
         if (right_pos < size())
            result = result && compare_(item.second, items_[heap_[right_pos]].second);
      }
   }

   return result;
}

// --------------------------------------------------------------------------------------------------------------------------------
