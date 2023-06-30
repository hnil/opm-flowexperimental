// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Contains the classes required to extend the black-oil model by energy.
 */
#ifndef EWOMS_BLACK_OIL_ENERGY_MODULE_TPFA_HH
#define EWOMS_BLACK_OIL_ENERGY_MODULE_TPFA_HH

#include "blackoilproperties.hh"
#include <opm/models/io/vtkblackoilenergymodule.hh>
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/common/Tabulated1DFunction.hpp>

#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

#include <string>

namespace Opm {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by energy.
 */
template <class TypeTag, bool enableEnergyV = getPropValue<TypeTag, Properties::EnableEnergy>()>
class BlackOilEnergyModule
