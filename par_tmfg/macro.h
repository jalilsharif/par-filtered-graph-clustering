#pragma once


#ifndef SIZE_T_MAX
#define SIZE_T_MAX std::numeric_limits<size_t>::max()
#endif
#define DBHT_NONE std::numeric_limits<std::size_t>::max()
#define NO_REACH false
#define NO_ATTACH_MIN std::numeric_limits<T>::min()
#define NO_ATTACH_MAX std::numeric_limits<T>::max()
#define NO_NEIGH SIZE_T_MAX