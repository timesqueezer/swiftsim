/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#include "adiabatic_index.h"
#include "equation_of_state.h"
#include "hydro_gradients.h"
#include "hydro_slope_limiters.h"
#include "io_properties.h"
#include "riemann.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
void hydro_read_particles(struct part* parts, struct io_props* list,
                          int* num_fields) {

  *num_fields = 8;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, parts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                parts, conserved.mass);
  list[3] = io_make_input_field("SmoothingLength", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, h);
  list[4] = io_make_input_field("InternalEnergy", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                conserved.energy);
  list[5] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_input_field("Accelerations", FLOAT, 3, OPTIONAL,
                                UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] = io_make_input_field("Density", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_DENSITY, parts, primitives.rho);
}

/**
 * @brief Get the internal energy of a particle
 *
 * @param e #engine.
 * @param p Particle.
 * @return Internal energy of the particle
 */
float convert_u(struct engine* e, struct part* p) {
  if (p->primitives.rho > 0.) {
    return gas_internal_energy_from_pressure(p->primitives.rho,
                                             p->primitives.P);
  } else {
    return 0.;
  }
}

/**
 * @brief Get the entropic function of a particle
 *
 * @param e #engine.
 * @param p Particle.
 * @return Entropic function of the particle
 */
float convert_A(struct engine* e, struct part* p) {
  if (p->primitives.rho > 0.) {
    return gas_entropy_from_pressure(p->primitives.rho, p->primitives.P);
  } else {
    return 0.;
  }
}

/**
 * @brief Get the total energy of a particle
 *
 * @param e #engine.
 * @param p Particle.
 * @return Total energy of the particle
 */
float convert_Etot(struct engine* e, struct part* p) {
#ifdef SHADOWFAX_TOTAL_ENERGY
  return p->conserved.energy;
#else
  if (p->conserved.mass > 0.) {
    float momentum2;

    momentum2 = p->conserved.momentum[0] * p->conserved.momentum[0] +
                p->conserved.momentum[1] * p->conserved.momentum[1] +
                p->conserved.momentum[2] * p->conserved.momentum[2];

    return p->conserved.energy + 0.5f * momentum2 / p->conserved.mass;
  } else {
    return 0.;
  }
#endif
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
void hydro_write_particles(struct part* parts, struct io_props* list,
                           int* num_fields) {

  *num_fields = 13;

  /* List what we want to write */
  list[0] = io_make_output_field("Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH,
                                 parts, x);
  list[1] = io_make_output_field("Velocities", FLOAT, 3, UNIT_CONV_SPEED, parts,
                                 primitives.v);
  list[2] = io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, parts,
                                 conserved.mass);
  list[3] = io_make_output_field("SmoothingLength", FLOAT, 1, UNIT_CONV_LENGTH,
                                 parts, h);
  list[4] = io_make_output_field_convert_part("InternalEnergy", FLOAT, 1,
                                              UNIT_CONV_ENERGY_PER_UNIT_MASS,
                                              parts, primitives.P, convert_u);
  list[5] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                 UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_output_field("Acceleration", FLOAT, 3,
                                 UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] = io_make_output_field("Density", FLOAT, 1, UNIT_CONV_DENSITY, parts,
                                 primitives.rho);
  list[8] = io_make_output_field("Volume", FLOAT, 1, UNIT_CONV_VOLUME, parts,
                                 cell.volume);
  list[9] = io_make_output_field("GradDensity", FLOAT, 3, UNIT_CONV_DENSITY,
                                 parts, primitives.gradients.rho);
  list[10] = io_make_output_field_convert_part(
      "Entropy", FLOAT, 1, UNIT_CONV_ENTROPY, parts, primitives.P, convert_A);
  list[11] = io_make_output_field("Pressure", FLOAT, 1, UNIT_CONV_PRESSURE,
                                  parts, primitives.P);
  list[12] =
      io_make_output_field_convert_part("TotEnergy", FLOAT, 1, UNIT_CONV_ENERGY,
                                        parts, conserved.energy, convert_Etot);
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
void writeSPHflavour(hid_t h_grpsph) {
  /* Gradient information */
  io_write_attribute_s(h_grpsph, "Gradient reconstruction model",
                       HYDRO_GRADIENT_IMPLEMENTATION);

  /* Slope limiter information */
  io_write_attribute_s(h_grpsph, "Cell wide slope limiter model",
                       HYDRO_SLOPE_LIMITER_CELL_IMPLEMENTATION);
  io_write_attribute_s(h_grpsph, "Piecewise slope limiter model",
                       HYDRO_SLOPE_LIMITER_FACE_IMPLEMENTATION);

  /* Riemann solver information */
  io_write_attribute_s(h_grpsph, "Riemann solver type",
                       RIEMANN_SOLVER_IMPLEMENTATION);
}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
int writeEntropyFlag() { return 0; }