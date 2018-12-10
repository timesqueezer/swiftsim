/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Josh Borrow (joshua.borrow@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_VISCOSITY_H
#define SWIFT_VISCOSITY_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"
#include "hydro_properties.h"
#include "kernel_hydro.h"
#include "part.h"

/* Choose the right functions */
#if defined(CONSTANT_VISCOSITY)
#include "./viscosity/Constant/viscosity.h"
#define VISCOSITY_IMPLEMENTATION "Basic constant artificial viscosity"
#endif


#endif /* SWIFT_VISCOSITY_H */