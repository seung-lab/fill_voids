/*
 * This file is part of fill_voids.
 * 
 * fill_voids is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * fill_voids is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with fill_voids.  If not, see <https://www.gnu.org/licenses/>.
 *
 * 
 * Author: William Silversmith
 * Affiliation: Seung Lab, Princeton University
 * Date: December 2019
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <stack>
#include <string>

#include "libdivide.h"

#ifndef FILLVOIDS_HPP
#define FILLVOIDS_HPP

#define BACKGROUND 0
#define VISITED_BACKGROUND 1
#define FOREGROUND 2

namespace fill_voids {

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
      yminus = yminus || (visited[cur-sx] == FOREGROUND);
    }
    else if (yminus) {
      stack.push( cur - sx );
      yminus = false;
    }
  }

  if (y < sy - 1) {
    if (visited[cur+sx]) {
      yplus = yplus || (visited[cur+sx] == FOREGROUND);
    }
    else if (yplus) {
      stack.push( cur + sx );
      yplus = false;
    }
  }

  if (z > 0) {
    if (visited[cur-sxyv]) {
      zminus = zminus || (visited[cur-sxyv] == FOREGROUND);
    }
    else if (zminus) {
      stack.push( cur - sxyv );
      zminus = false;
    }
  }

  if (z < sz - 1) {
    if (visited[cur+sxyv]) {
      zplus = zplus || (visited[cur+sxyv] == FOREGROUND);
    }
    else if (zplus) {
      stack.push( cur + sxyv );
      zplus = false;
    }
  }
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
void _binary_fill_holes(
  T* labels, 
  const size_t sx, const size_t sy, const size_t sz
) {

  const size_t sxy = sx * sy;
  const size_t voxels = sx * sy * sz;

  if (voxels == 0) {
    return;
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
      labels[cur] = VISITED_BACKGROUND;
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
      labels[cur] = VISITED_BACKGROUND;
      add_neighbors<T>(
        labels, stack,
        sx, sy, sz, 
        cur, y, z,
        yplus, yminus, zplus, zminus
      );
    }    
  }

  for (size_t i = 0; i < voxels; i++) {
    labels[i] = static_cast<T>(labels[i] != VISITED_BACKGROUND);
  }
}


};

#undef BACKGROUND
#undef VISITED_BACKGROUND
#undef FOREGROUND

#endif