/*******************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - header file
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    May 01 2010  v.0.1    first release
              May 10 2014  v.0.2    solution, option and info structures updated, new prototypes added
              Jun 28 2019  v.0.3.0  new and more efficient organization of data structures and functions
              Mar 21 2020  v.0.3.1  adding triplet structure and new function prototypes
              May 19 2020  v.0.3.2  reorganization of OPTION structure, new function prototypes
              Apr 13 2022  v.0.3.2  patch
********************************************************************************************************/
#ifndef _BP_
#define _BP_
#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>
#include <array>
#include <map>
#include <chrono>
#include <memory>
#include <typeinfo>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <signal.h>
#include <sys/time.h>
#endif

extern "C"
{
#include "bp.h"
}