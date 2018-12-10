/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Josh Borrow (joshua.borrow@durham.ac.uk)
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
#ifndef SWIFT_CONSTANT_VISCOSITY_PART_H
#define SWIFT_CONSTANT_VISCOSITY_PART_H

/**
 * @file Constant/viscosity_part.h
 * @brief Struct for viscosity related quantities in the basic (constant) version
 *        of artificial viscosities.
 * 
 * These terms are the ones expressed in the first AV paper by Monaghan,
 * Journal of Computational Physics Volume 136, Issue 2, 15 September 1997,
 * Pages 298-307.
 * 
 * We use the same numerical coefficients as the Gadget-2 code.
 */

/* Viscosity information for given particle */
struct viscosity_part {

  /* TODO: Union this */

  struct {

  } density;

  struct {

  } force;

} SWIFT_STRUCT_ALIGN;

/* Viscosity information for hydro_properties */
struct viscosity_properties {

  /* Viscosity coefficient (constant) */
  float alpha;

};

#endif /* SWIFT_CONSTANT_VISCOSITY_PART_H */