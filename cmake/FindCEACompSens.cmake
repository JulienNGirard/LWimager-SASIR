# $Id: FindCEA_COMP_SENS.cmake 13814 2009-08-20 11:55:06Z loose $
#
# Copyright (C) 2008-2009
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, seg@astron.nl
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Try to find CEA_COMP_SENS.
#
# Variables used by this module:
#  CEA_COMP_SENS_ROOT_DIR     - CEA_COMP_SENS root directory
#
# Variables defined by this module:
#  CEA_COMP_SENS_FOUND        - system has CEA_COMP_SENS
#  CEA_COMP_SENS_INCLUDE_DIR  - the CEA_COMP_SENS include directory (cached)
#  CEA_COMP_SENS_INCLUDE_DIRS - the CEA_COMP_SENS include directories
#                         (identical to CEA_COMP_SENS_INCLUDE_DIR)
#  CEA_COMP_SENS_LIBRARY      - the CEA_COMP_SENS library (cached)
#  CEA_COMP_SENS_LIBRARIES    - the CEA_COMP_SENS libraries
#                         (identical to CEA_COMP_SENS_LIBRARY)

if(NOT CEA_COMP_SENS_FOUND)

  find_path(CEA_COMP_SENS_INCLUDE_DIR CEA_comp_sens.h
    HINTS ${CEA_COMP_SENS_ROOT_DIR} PATH_SUFFIXES include)
  find_library(CEA_COMP_SENS_LIBRARY CEA_comp_sens
    HINTS ${CEA_COMP_SENS_ROOT_DIR} PATH_SUFFIXES lib)
  find_library(M_LIBRARY m)
  mark_as_advanced(CEA_COMP_SENS_INCLUDE_DIR CEA_COMP_SENS_LIBRARY M_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(CEA_COMP_SENS DEFAULT_MSG
    CEA_COMP_SENS_LIBRARY M_LIBRARY CEA_COMP_SENS_INCLUDE_DIR)

  set(CEA_COMP_SENS_INCLUDE_DIRS ${CEA_COMP_SENS_INCLUDE_DIR})
  set(CEA_COMP_SENS_LIBRARIES ${CEA_COMP_SENS_LIBRARY} ${M_LIBRARY})

endif(NOT CEA_COMP_SENS_FOUND)
