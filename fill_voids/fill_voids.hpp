/*
 * This file is part of fill_voids.
 * 
 * fill_voids is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * fill_voids is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 * 
 * You should have received a copy of the Lesser GNU General Public License
 * along with fill_voids.  If not, see <https://www.gnu.org/licenses/>.
 *
 * 
 * Author: William Silversmith
 * Affiliation: Seung Lab, Princeton University
 * Date: December 2019 - Februrary 2020
 */

#ifndef FILLVOIDS_HPP
#define FILLVOIDS_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <memory>
#include <vector>
#include <stack>
#include <string>

#include "libdivide.h"

namespace fill_voids {

enum Label {
  BACKGROUND = 0,
  VISITED_BACKGROUND = 1,
  FOREGROUND = 2
};

template <typename T>
class DisjointSet {
public:
  T *ids;
  size_t length;

  DisjointSet () {
    length = 65536; // 2^16, some "reasonable" starting size
    ids = new T[length]();
  }

  DisjointSet (size_t len) {
    length = len;
    ids = new T[length]();
  }

  DisjointSet (const DisjointSet &cpy) {
    length = cpy.length;
    ids = new T[length]();

    for (int i = 0; i < length; i++) {
      ids[i] = cpy.ids[i];
    }
  }

  ~DisjointSet () {
    delete []ids;
  }

  T root (T n) {
    T i = ids[n];
    while (i != ids[i]) {
      ids[i] = ids[ids[i]]; // path compression
      i = ids[i];
    }

    return i;
  }

  bool find (T p, T q) {
    return root(p) == root(q);
  }

  void add(T p) {
    if (p >= length) {
      printf("Connected Components Error: Label %lli cannot be mapped to union-find array of length %lu.\n", static_cast<long long int>(p), length);
      throw std::runtime_error("maximum length exception");
    }

    if (ids[p] == 0) {
      ids[p] = p;
    }
  }

  void unify (T p, T q) {
    if (p == q) {
      return;
    }

    T i = root(p);
    T j = root(q);

    if (i == 0) {
      add(p);
      i = p;
    }

    if (j == 0) {
      add(q);
      j = q;
    }

    if (i < j) {
      ids[j] = i;
    }
    else {
      ids[i] = j;
    }
  }

  void clear() {
    for (int64_t i = 0; i < length; i++) {
      ids[i] = i;
    }
  }

  void print() {
    for (int i = 0; i < length; i++) {
      printf("%d, ", ids[i]);
    }
    printf("\n");
  }

  // would be easy to write remove. 
  // Will be O(n).
};

uint8_t online_relabel(
  uint8_t* out_labels, 
  const int64_t sx, const int64_t sy, const int64_t sz,
  const int64_t start_loc, const int64_t end_loc,
  DisjointSet<uint8_t>& equivalences,
  const std::vector<uint8_t>& current_active_labels,
  const std::vector<uint8_t>& previous_active_labels
) {
  uint8_t label = 0;
  uint8_t renumber[256] = {};
  bool active[256] = {};

  for (const uint8_t i : current_active_labels) {
    active[i] = true;
  }

  for (const uint8_t i : previous_active_labels) {
    active[i] = true;
  }

  uint8_t next_cannonical_label = 2;

  for (int64_t i = 1; i < 256; i++) {
    label = equivalences.root(i);
    if (label <= 1) {
      renumber[i] = label;
    }
    else if (label & !active[label]) {
      renumber[i] = 1;
      equivalences.ids[i] = 1;
    }
    else {
      renumber[i] = next_cannonical_label++;
    }
  }

  // Raster Scan 2: Write final labels based on equivalences
  for (int64_t i = start_loc; i < end_loc; i++) {
    out_labels[i] = renumber[out_labels[i]];
  }

  return next_cannonical_label;
}

uint8_t final_relabel(
  uint8_t* out_labels, 
  const int64_t sx, const int64_t sy, const int64_t sz,
  const int64_t start_loc, const int64_t end_loc,
  DisjointSet<uint8_t>& equivalences,
  const std::vector<uint8_t>& current_active_labels,
  const std::vector<uint8_t>& previous_active_labels
) {
  uint8_t label = 0;
  uint8_t renumber[256] = {};
  bool active[256] = {};

  for (const uint8_t i : current_active_labels) {
    active[i] = true;
  }

  for (const uint8_t i : previous_active_labels) {
    active[i] = true;
  }

  uint8_t next_cannonical_label = 2;

  for (int64_t i = 1; i < 256; i++) {
    label = equivalences.root(i);
    if (label <= 1) {
      renumber[i] = label;
    }
    else if (label & !active[label]) {
      renumber[i] = 1;
      equivalences.ids[i] = 1;
    }
    else {
      renumber[i] = next_cannonical_label++;
    }
  }

  // Raster Scan 2: Write final labels based on equivalences
  for (int64_t i = start_loc; i < end_loc; i++) {
    out_labels[i] = renumber[out_labels[i]];
  }

  return next_cannonical_label;
}


template <typename T>
int64_t binary_fill_holes3d_ccl(
    T* in_labels, 
    const int64_t sx, const int64_t sy, const int64_t sz,
    uint8_t* out_labels
) {

  const int64_t sxy = sx * sy;
  const int64_t voxels = sxy * sz;

  DisjointSet<uint8_t> equivalences(256);

  equivalences.add(0);
  equivalences.add(1);

  /*
    Layout of forward pass mask (which faces backwards). 
    N is the current location.

    z = -1     z = 0
    A B C      J K L   y = -1 
    D E F      M N     y =  0
    G H I              y = +1
   -1 0 +1    -1 0   <-- x axis
  */

  // Z - 1
  const int64_t B = -sx - sxy;
  const int64_t E = -sxy;
  const int64_t D = -1 - sxy;

  // Current Z
  const int64_t K = -sx;
  const int64_t M = -1;
  const int64_t J = -1 - sx;
  // N = 0;

  int64_t loc = 0;
  uint8_t next_label = 1;
  int64_t relabel_from = 0;
  bool touching_border = true;

  std::vector<uint8_t> previous_active_labels;
  std::vector<uint8_t> current_active_labels;
  current_active_labels.reserve(256);
  previous_active_labels.reserve(256);

  int64_t original_foreground_count = 0;

  equivalences.print();

  // Raster Scan 1: Set temporary labels and 
  // record equivalences in a disjoint set.

  for (int64_t z = 0; z < sz; z++) {
    std::swap(current_active_labels, previous_active_labels);
    current_active_labels.clear();

    for (int64_t y = 0; y < sy; y++) {
      // printf("%d, %d | %d\n", y, z, next_label);  
      const int64_t row_start = sx * (y + sy * z);
      loc = row_start;
      T cur = in_labels[row_start];

      original_foreground_count += (cur != 0);

      if (cur != 0) {
        if (y > 0 && cur == in_labels[loc + K]) {
          out_labels[loc] = out_labels[loc + K];
          if (z > 0 && cur == in_labels[loc + E]) {
            equivalences.unify(out_labels[loc], out_labels[loc + E]);
          }
        }
        else if (z > 0 && cur == in_labels[loc + E]) {
          out_labels[loc] = out_labels[loc + E];
        }
        else if (next_label == 255) {
          // equivalences.print();
          next_label = online_relabel(out_labels, sx, sy, sz, relabel_from, loc, equivalences, current_active_labels, previous_active_labels);
          out_labels[loc] = next_label;
          equivalences.clear();
          equivalences.add(next_label);
          current_active_labels.push_back(next_label);
          relabel_from = loc;
        }
        else {
          next_label++;
          out_labels[loc] = next_label;
          equivalences.add(out_labels[loc]);
          current_active_labels.push_back(next_label);
        }
      }

      touching_border = cur == 0;

      for (int64_t x = 1; x < sx; x++) {
        loc = row_start + x;
        cur = in_labels[loc];

        original_foreground_count += (cur > 0);

        if (cur) {
          touching_border = false;
          out_labels[loc] = 1;
          continue;
        }
        else if (touching_border) {
          continue;
        }

        if (x > 0 && cur == in_labels[loc + M]) {
          out_labels[loc] = out_labels[loc + M];

          if (y > 0 && cur == in_labels[loc + K] && cur != in_labels[loc + J]) {
            equivalences.unify(out_labels[loc], out_labels[loc + K]); 
            if (z > 0 && cur == in_labels[loc + E]) {
              if (cur != in_labels[loc + D] && cur != in_labels[loc + B]) {
                equivalences.unify(out_labels[loc], out_labels[loc + E]);
              }
            }
          }
          else if (z > 0 && cur == in_labels[loc + E] && cur != in_labels[loc + D]) {
            equivalences.unify(out_labels[loc], out_labels[loc + E]); 
          }
        }
        else if (y > 0 && cur == in_labels[loc + K]) {
          out_labels[loc] = out_labels[loc + K];

          if (z > 0 && cur == in_labels[loc + E] && cur != in_labels[loc + B]) {
            equivalences.unify(out_labels[loc], out_labels[loc + E]); 
          }
        }
        else if (z > 0 && cur == in_labels[loc + E]) {
          out_labels[loc] = out_labels[loc + E];
        }
        else if (next_label == 255) {
          // equivalences.print();
          next_label = online_relabel(out_labels, sx, sy, sz, relabel_from, loc, equivalences, current_active_labels, previous_active_labels);
          // printf("next label: %d\n", next_label);
          out_labels[loc] = next_label;
          // printf("wow\n");
          equivalences.clear();
          // printf("cleared\n");
          equivalences.add(next_label);
          current_active_labels.push_back(next_label);
          // printf("add\n");
          relabel_from = loc;
          // printf("relabel_from\n");
        }
        else {
          next_label++;
          out_labels[loc] = next_label;
          equivalences.add(out_labels[loc]);
          current_active_labels.push_back(next_label);
        }
      }
    }
  }

  // printf("nl %d\n", next_label);

  for (int64_t z = 0; z < sz; z++) {
    for (int64_t y = 0; y < sy; y++) {
      loc = sx * (y + sy * z);
      if (in_labels[loc] == 0) {
        equivalences.unify(0, out_labels[loc]);
      }
      loc = (sx-1) + sx * (y + sy * z);
      if (in_labels[loc] == 0) {
        equivalences.unify(0, out_labels[loc]);
      }
    }
  }

  for (int64_t z = 0; z < sz; z++) {
    for (int64_t x = 0; x < sx; x++) {
      loc = x + sxy * z;
      if (in_labels[loc] == 0) {
        equivalences.unify(0, out_labels[loc]);
      }
      loc = x + sx * (sy-1) + sxy * z;
      if (in_labels[loc] == 0) {
        equivalences.unify(0, out_labels[loc]);
      }
    }
  }
  for (int64_t y = 0; y < sy; y++) {
    for (int64_t x = 0; x < sx; x++) {
      loc = x + sx * y;
      if (in_labels[loc] == 0) {
        equivalences.unify(0, out_labels[loc]);
      }
      loc = x + sx * y + sxy * (sz-1);
      if (in_labels[loc] == 0) {
        equivalences.unify(0, out_labels[loc]);
      }
    }
  }

  // printf("relabel_from: %d\n", relabel_from);

  // printf("equivalences:\n");
  // equivalences.print();

  for (int i = 2; i < 256; i++) {
    equivalences.ids[i] = (equivalences.root(i) > 0);
  }

  // printf("equivalences:\n");
  // equivalences.print();

  online_relabel(out_labels, sx, sy, sz, relabel_from, voxels, equivalences, current_active_labels, previous_active_labels);

  int64_t final_foreground_count = 0;
  for (int64_t i = 0; i < voxels; i++) {
    final_foreground_count += out_labels[i] > 0;
  }

  const int64_t num_filled = final_foreground_count - original_foreground_count;

  return num_filled;
}

template <typename T>
inline void push_stack(
  T* labels, const size_t loc,
  std::stack<size_t> &stack, bool &placed
) {
  if (labels[loc] == 0) {
    if (!placed) {
      stack.push(loc);
    }
    placed = true;
  }
  else {
    placed = false;
  }  
}

template <typename T>
inline void add_neighbors(
  T* visited, std::stack<size_t> &stack,
  const size_t sx, const size_t sy,
  const size_t cur, const size_t y,
  bool &yplus, bool &yminus
) {
  // Only add a seed point if we've just 
  // started OR have just passed a foreground
  // voxel.

  if (y > 0) {
    if (visited[cur-sx]) {
      yminus = yminus || (visited[cur-sx] == Label::FOREGROUND);
    }
    else if (yminus) {
      stack.push( cur - sx );
      yminus = false;
    }
  }

  if (y < sy - 1) {
    if (visited[cur+sx]) {
      yplus = yplus || (visited[cur+sx] == Label::FOREGROUND);
    }
    else if (yplus) {
      stack.push( cur + sx );
      yplus = false;
    }
  }
}

template <typename T>
inline void add_neighbors(
  T* visited, std::stack<size_t> &stack,
  const size_t sx, const size_t sy, const size_t sz, 
  const size_t cur, const size_t y, const size_t z,
  bool &yplus, bool &yminus, bool &zplus, bool &zminus
) {
  const size_t sxyv = sx * sy;

  // Only add a seed point if we've just 
  // started OR have just passed a foreground
  // voxel.

  if (y > 0) {
    if (visited[cur-sx]) {
      yminus = yminus || (visited[cur-sx] == Label::FOREGROUND);
    }
    else if (yminus) {
      stack.push( cur - sx );
      yminus = false;
    }
  }

  if (y < sy - 1) {
    if (visited[cur+sx]) {
      yplus = yplus || (visited[cur+sx] == Label::FOREGROUND);
    }
    else if (yplus) {
      stack.push( cur + sx );
      yplus = false;
    }
  }

  if (z > 0) {
    if (visited[cur-sxyv]) {
      zminus = zminus || (visited[cur-sxyv] == Label::FOREGROUND);
    }
    else if (zminus) {
      stack.push( cur - sxyv );
      zminus = false;
    }
  }

  if (z < sz - 1) {
    if (visited[cur+sxyv]) {
      zplus = zplus || (visited[cur+sxyv] == Label::FOREGROUND);
    }
    else if (zplus) {
      stack.push( cur + sxyv );
      zplus = false;
    }
  }
}

/* The idea here is to scan the four sides
 * of the image and insert an exploration
 * point into the stack whenever an exterior  
 * void (defined as touching the edge of the 
 * image) is first encountered.
 */
template <typename T>
void initialize_stack(
    T* labels, 
    const size_t sx, const size_t sy,
    std::stack<size_t> &stack
  ) {

  bool placed_front = false;
  bool placed_back = false;

  size_t loc;
  for (size_t x = 0; x < sx; x++) {
    loc = x;
    push_stack<T>(labels, loc, stack, placed_front);
    
    loc = x + sx * (sy - 1);
    push_stack<T>(labels, loc, stack, placed_back);
  }

  placed_front = false;
  placed_back = false;

  for (size_t y = 0; y < sy; y++) {
    loc = sx * y;
    push_stack<T>(labels, loc, stack, placed_front);
    
    loc = (sx - 1) + sx * y;
    push_stack<T>(labels, loc, stack, placed_back);
  }
}

/* The idea here is to scan the six faces
 * of the image and insert an exploration
 * point into the stack whenever an exterior  
 * void (defined as touching the edge of the 
 * image) is first encountered.
 *
 * This is a lower memory version of the previous
 * logic, which copied the entire image into a buffer 
 * with a 1 pixel black border and explored from <0,0,0> 
 * which exploits the knowledge that that border touches
 * everything and will find those exterior voids automatically.
 */
template <typename T>
void initialize_stack(
    T* labels, 
    const size_t sx, const size_t sy, const size_t sz,
    std::stack<size_t> &stack
  ) {
  const size_t sxy = sx * sy;

  bool placed_front = false;
  bool placed_back = false;

  size_t loc;
  for (size_t y = 0; y < sy; y++) {
    for (size_t x = 0; x < sx; x++) {
      loc = x + sx * y;
      push_stack<T>(labels, loc, stack, placed_front);
      
      loc = x + sx * y + sxy * (sz - 1);
      push_stack<T>(labels, loc, stack, placed_back);
    }
  }

  placed_front = false;
  placed_back = false;

  for (size_t z = 0; z < sz; z++) {
    for (size_t x = 0; x < sx; x++) {
      loc = x + sxy * z;
      push_stack<T>(labels, loc, stack, placed_front);
      
      loc = x + sx * (sy - 1) + sxy * z;
      push_stack<T>(labels, loc, stack, placed_back);
    }
  }

  placed_front = false;
  placed_back = false;

  for (size_t z = 0; z < sz; z++) {
    for (size_t y = 0; y < sy; y++) {
      loc = sx * y + sxy * z;
      push_stack<T>(labels, loc, stack, placed_front);

      loc = (sx - 1) + sx * y + sxy * z;
      push_stack<T>(labels, loc, stack, placed_back); 
    }
  }
}

template <typename T>
size_t binary_fill_holes2d(
  T* labels, 
  const size_t sx, const size_t sy
) {
  
  const size_t voxels = sx * sy;

  if (voxels == 0) {
    return 0;
  }

  // mark all foreground as 2 (FOREGROUND) 
  // so we can mark visited as 1 (VISITED_BACKGROUND) 
  // without overwriting foreground as we want foreground 
  // to be 2 and voids to be 0 (BACKGROUND)
  for (size_t i = 0; i < voxels; i++) {
    labels[i] = static_cast<T>(static_cast<uint8_t>(labels[i] != 0) * 2);
  }

  const libdivide::divider<size_t> fast_sx(sx); 

  std::stack<size_t> stack; 
  initialize_stack(labels, sx, sy, stack);

  while (!stack.empty()) {
    size_t loc = stack.top();
    stack.pop();

    if (labels[loc]) {
      continue;
    }

    size_t y = loc / fast_sx;
    size_t startx = y * sx;

    bool yplus = true;
    bool yminus = true;

    for (size_t cur = loc; cur < startx + sx; cur++) {
      if (labels[cur]) {
        break;
      }
      labels[cur] = Label::VISITED_BACKGROUND;
      add_neighbors<T>(
        labels, stack,
        sx, sy, 
        cur, y,
        yplus, yminus
      );
    }

    yplus = true;
    yminus = true;

    // avoid integer underflow
    for (int64_t cur = static_cast<int64_t>(loc) - 1; cur >= static_cast<int64_t>(startx); cur--) {
      if (labels[cur]) {
        break;
      }
      labels[cur] = Label::VISITED_BACKGROUND;
      add_neighbors<T>(
        labels, stack,
        sx, sy,
        cur, y,
        yplus, yminus
      );
    }    
  }

  size_t num_filled = 0;
  for (size_t i = 0; i < voxels; i++) {
    num_filled += static_cast<size_t>(labels[i] == Label::BACKGROUND);
    labels[i] = static_cast<T>(labels[i] != Label::VISITED_BACKGROUND);
  }

  return num_filled;
}

template <typename T>
size_t binary_fill_holes3d(
  T* labels, 
  const size_t sx, const size_t sy, const size_t sz
) {

  const size_t sxy = sx * sy;
  const size_t voxels = sx * sy * sz;

  if (voxels == 0) {
    return 0;
  }

  // mark all foreground as 2 (FOREGROUND) 
  // so we can mark visited as 1 (VISITED_BACKGROUND) 
  // without overwriting foreground as we want foreground 
  // to be 2 and voids to be 0 (BACKGROUND)
  for (size_t i = 0; i < voxels; i++) {
    labels[i] = static_cast<T>(static_cast<uint8_t>(labels[i] != 0) * 2);
  }

  const libdivide::divider<size_t> fast_sx(sx); 
  const libdivide::divider<size_t> fast_sxy(sxy); 

  std::stack<size_t> stack; 
  initialize_stack(labels, sx, sy, sz, stack);

  while (!stack.empty()) {
    size_t loc = stack.top();
    stack.pop();

    if (labels[loc]) {
      continue;
    }

    size_t z = loc / fast_sxy;
    size_t y = (loc - (z * sxy)) / fast_sx;
    size_t startx = y * sx + z * sxy;

    bool yplus = true;
    bool yminus = true;
    bool zplus = true;
    bool zminus = true;

    for (size_t cur = loc; cur < startx + sx; cur++) {
      if (labels[cur]) {
        break;
      }
      labels[cur] = Label::VISITED_BACKGROUND;
      add_neighbors<T>(
        labels, stack,
        sx, sy, sz, 
        cur, y, z,
        yplus, yminus, zplus, zminus
      );
    }

    yplus = true;
    yminus = true;
    zplus = true;
    zminus = true;

    // avoid integer underflow
    for (int64_t cur = static_cast<int64_t>(loc) - 1; cur >= static_cast<int64_t>(startx); cur--) {
      if (labels[cur]) {
        break;
      }
      labels[cur] = Label::VISITED_BACKGROUND;
      add_neighbors<T>(
        labels, stack,
        sx, sy, sz, 
        cur, y, z,
        yplus, yminus, zplus, zminus
      );
    }    
  }

  size_t num_filled = 0;
  for (size_t i = 0; i < voxels; i++) {
    num_filled += static_cast<size_t>(labels[i] == Label::BACKGROUND);
    labels[i] = static_cast<T>(labels[i] != Label::VISITED_BACKGROUND);
  }

  return num_filled;
}

template <typename T>
size_t binary_fill_holes(
  T* labels, 
  const size_t sx, const size_t sy, const size_t sz
) {
  return binary_fill_holes3d<T>(labels, sx, sy, sz);
}

template <typename T>
size_t binary_fill_holes(
  T* labels, 
  const size_t sx, const size_t sy
) {
  return binary_fill_holes2d<T>(labels, sx, sy);
}

};

#endif