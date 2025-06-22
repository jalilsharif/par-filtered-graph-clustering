#pragma once

#include <atomic>
#include <cstring>
#include <iostream>
#include <cstdlib>

// Enable MCX16 for x86 if required
#define MCX16

// https://raw.githubusercontent.com/cmuparlay/pbbsbench/37df3e14c7f3d738500e06840874a1505944598d/common/atomics.h

namespace pbbs {

  template <typename ET>
  inline bool atomic_compare_and_swap(ET* a, ET oldval, ET newval) {
#ifdef __aarch64__ // ARM 64-bit architecture
    if (sizeof(ET) == 1) {
      uint8_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 4) {
      uint32_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 8) {
      uint64_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 16) {
      __int128 r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __atomic_compare_exchange(reinterpret_cast<volatile __int128*>(a),
                                       &r_oval, &r_nval, false, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST);
    } else {
      std::cerr << "Unsupported CAS size for ARM: " << sizeof(ET) << " bytes" << std::endl;
      std::exit(EXIT_FAILURE);
    }
#elif defined(MCX16) // x86 with MCX16 support
    if (sizeof(ET) == 1) {
      uint8_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 4) {
      uint32_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 8) {
      uint64_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 16) {
      __int128 r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap_16(reinterpret_cast<volatile __int128*>(a), r_oval, r_nval);
    } else {
      std::cerr << "Unsupported CAS size for x86: " << sizeof(ET) << " bytes" << std::endl;
      std::exit(EXIT_FAILURE);
    }
#else
    std::cerr << "Unsupported platform for atomic_compare_and_swap" << std::endl;
    std::exit(EXIT_FAILURE);
#endif
  }

  template <typename E, typename EV>
  inline E fetch_and_add(E *a, EV b) {
    E newV, oldV;
    do {
      oldV = *a;
      newV = oldV + b;
    } while (!atomic_compare_and_swap(a, oldV, newV));
    return oldV;
  }

  template <typename E, typename EV>
  EV write_add(E *a, EV b) {
    E newV, oldV;
    do {
      oldV = *a;
      newV = oldV + b;
    } while (!atomic_compare_and_swap(a, oldV, newV));
    return newV;
  }

  template <typename ET, typename F>
  inline bool write_min(ET *a, ET b, F less) {
    ET c;
    bool r = 0;
    do c = *a;
    while (less(b, c) && !(r = atomic_compare_and_swap(a, c, b)));
    return r;
  }

  template <typename ET, typename F>
  inline bool write_max(ET *a, ET b, F less) {
    ET c;
    bool r = 0;
    do c = *a;
    while (less(c, b) && !(r = atomic_compare_and_swap(a, c, b)));
    return r;
  }
}
