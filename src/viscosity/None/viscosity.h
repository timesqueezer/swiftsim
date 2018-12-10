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
#ifndef SWIFT_NONE_VISCOSITY_H
#define SWIFT_NONE_VISCOSITY_H

/**
 * @file None/viscosity.h
 * @brief Empty viscosity interaction functions (used for schemes that do not
 *        require an artificial viscosity, such as GIZMO.)
 * 
 * These terms are the ones expressed in the first AV paper by Monaghan,
 * Journal of Computational Physics Volume 136, Issue 2, 15 September 1997,
 * Pages 298-307.
 * 
 * We use the same numerical coefficients as the Gadget-2 code.
 */

#include "hydro.h"
#include "hydro_properties.h"

/**
 * @brief Computes the hydro time-step of a given particle according to the
 *        viscosity scheme. This can be empty. If it is not empty, then the
 *        minimum of this and the result from hydro_compute_timestep will
 *        be used.
 * 
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 * @param hydro_properties The constants used in the scheme
 * @param cosmo The cosmological model.
 */
float viscosity_compute_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct hydro_props *restrict hydro_properties,
    const struct cosmology *restrict cosmo) {
  return FLT_MAX;    
}


/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various density tasks
 * 
 * @param p The particle to act upon
 * @param hs #hydro_space containing hydro specific space information.
 */
void viscosity_init_part(
    struct part *restrict p, const struct hydro_space *hs) {}


/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 * 
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
void viscosity_end_density(
    struct part *restrict p, const struct cosmology *cosmo) {}


/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 * 
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 */
void viscosity_part_has_no_neighbours(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo) {}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * This function is called in the ghost task to convert some quantities coming
 * from the density loop over neighbours into quantities ready to be used in the
 * force loop over neighbours. Quantities are typically read from the density
 * sub-structure and written to the force sub-structure.
 * Examples of calculations done here include the calculation of viscosity term
 * constants, thermal conduction terms, hydro conversions, etc.
 * 
 * This function is actually called within hydro_prepare_force to enable the
 * setting and getting of variables that are in the unions.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 * @param hydro_props Hydrodynamic properties.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 */
void viscosity_prepare_force(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const float dt_alpha) {}


/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the forces and accelerations by the appropiate constants
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
void viscosity_end_force(
    struct part *restrict p, const struct cosmology *cosmo) {}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
void viscosity_first_init_part(
    struct part *restrict p, struct xpart *restrict xp) {}

#endif /* SWIFT_NONE_VISCOSITY_H */