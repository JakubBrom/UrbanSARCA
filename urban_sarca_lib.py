#!/usr/bin/env python
# -*- coding: cp1250 -*- 

# ------------------------------------------------------------------------------
# Urban Green SARCA Library

# ------------------------------------------------------------------------------
#  Author: Jakub Brom
#  Date: 2020 - 10 - 07
#
#  Copyright (c)  Jakub Brom, 2020 - 2023.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------------

import numpy as np


def interceptFactor(LAI, precip, fresh_biomass, k=1.0, S=0.2):
	"""
	Interception factor for both dry and wet deposition of
	radionuclides.

	:param LAI: Leaf Area Index (unitless)
	:param k: Constant of radionuclide: I = 0.5, Sr and Ba = 2.0, Cs and \
	another radionuclides = 1.0
	:param precip: Precipitation amount (mm) for period of deposition \
	(ca 24 hours after radiation accident).
	:param fresh_biomass: Amount of fresh biomass :math:`(t.ha^{-1})`
	:param S: Mean thickness of water film at plant leaves (mm). \
	Default S =	0.2 mm

	:return: Interception Factor (rel.)
	"""

	try:
		IF = LAI * k * S * (1.0 - np.exp((-np.log(2.0))/(3.0 * S) * (
				precip + 0.0001))) / (precip + 0.0001)
		IF[IF > 1.0] = 1.0
		IF[fresh_biomass < 0.5] = 0.0
	except ArithmeticError:
		raise ArithmeticError("Interception factor has not been "
							  "calculated")

	return IF

def contBiomass(depo, IF):
	"""
	Radiaoctive contamination of biomass :math:`(Bq.m^{-2})`

	:param depo: Total radioactive deposition :math:`(Bq.m^{-2})`
	:param IF: Interception Factor (rel.)
	:return: Radioactive contamination of biomass :math:`(Bq.m^{-2})`
	"""

	try:
		cont_biomass = depo * IF
		cont_biomass[cont_biomass < 0.0] = 0.0
	except ArithmeticError:
		raise ArithmeticError("Vegetation biomass radioactive "
							  "contamination has not been calculated")

	return cont_biomass
	
def contSoil(depo, IF):
	"""
	Radiaoctive contamination of soil :math:`(Bq.m^{-2})`

	:param depo: Total radioactive deposition :math:`(Bq.m^{-2})`
	:param IF: Interception Factor (rel.)
	:return: Radioactive contamination of soil :math:`(Bq.m^{-2})`
	"""

	try:
		cont_soil = depo * (1 - IF)
		cont_soil[cont_soil < 0.0] = 0.0
	except ArithmeticError:
		raise ArithmeticError("Soil radioactive "
							  "contamination has not been calculated")

	return cont_soil

def contMass(cont_biomass, fresh_biomass):
	"""
	Calculation radioactive contamination of fresh vegetation mass \
	:math:`(Bq.kg^{-1})`

	:param cont_biomass: Radioactive deposition on biomass \
	:math:`(Bq.m^{-2})`
	:param fresh_biomass: Amount of fresh biomass of vegetation \
	:math:`(t.ha^{-1})`

	:return: Fresh vegetation mass radioactive contamination \
	:math:`(Bq.kg^{-1})`
	"""

	ignore_zero = np.seterr(all="ignore")

	try:
		cont_mass = cont_biomass/(fresh_biomass * 0.1)
		cont_mass[cont_mass < 0.0] = 0.0
	except ArithmeticError:
		raise ArithmeticError("Mass Radioactive "
							  "contamination of fresh biomass has not "
							  "been calculated")

	return cont_mass

def referLevel(depo, rl_dict):
	"""
	Mask of radioactive deposition for reference levels (categories) defined by tresholds.

	:param depo: Total radioactive deposition :math:`(Bq.m^{-2})`
	:param rl_dict: Dictionary with reference level tresholds :math:`(Bq.m^{-2})` and their keys (RL 1, RL 2, ...).
	:return: Mask of radioactive deposition for reference levels.
	"""

	try:
		# Sort dict according to values and extract RL values for categories
		rl_dict_sort = {k: v for k, v in sorted(rl_dict.items(), key=lambda item: item[1])}
		rl_keys = [int(i.split(' ')[1]) for i in list(rl_dict_sort.keys())]
		rl_values = list(rl_dict_sort.values())

		# Make ref_levels array
		reference_groups = np.where(depo < rl_values[0], 0.0, depo)
		for i in range(0,len(rl_keys)):
			reference_groups = np.where(depo >= rl_values[i], rl_keys[i], reference_groups)
	except ArithmeticError:
		raise ArithmeticError("Mask for reference levels has not been"
							  " calculated")

	return reference_groups

def hygLimit(cont_mass, hyg_lim=1000):
	"""
	Mask for hygienic limit of biomass radioactive contamination 
	exceeding.
	
	:param cont_mass: Fresh vegetation mass radioactive contamination \
	:math:`(Bq.kg^{-1})`
	:param hyg_lim: Hygienic limit of mass radioactive contamination \
	:math:`(Bq.kg^{-1})`. Default value is 1000 :math:`(Bq.kg^{-1})` \
	according to Czech law.
	:return: Mask for hygienic limit of biomass radioactive \
	contamination exceeding.
	"""

	try:
		hyg_lim_mask = np.where(cont_mass >= hyg_lim, 1.0, 0.0)
	except ArithmeticError:
		raise ArithmeticError("Mask for hygienic limit of mass "
							  "radioactive contamination has not been "
							  "calculated")

	return hyg_lim_mask