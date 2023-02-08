/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <marty/models/gthdm.h>
#include <marty/models/nmfv.h>
#include <marty/models/pmssm.h>
#include <marty/models/pmssm_lem.h>
#include <marty/models/qcd.h>
#include <marty/models/qed.h>
#include <marty/models/sm.h>

#include "CepGenAddOns/MartyWrapper/MartyModelFactory.h"

typedef mty::GTHDM_Model GTHDM;
REGISTER_MARTY_MODEL("gthdm", GTHDM)

//typedef mty::NMFV_Model NMFV;
//REGISTER_MARTY_MODEL("nmfv", NMFV)

typedef mty::PMSSM_Model PMSSM;
REGISTER_MARTY_MODEL("pmssm", PMSSM)

typedef mty::PMSSM_LEM PMSSM_LEM;
REGISTER_MARTY_MODEL("pmssm_lem", PMSSM_LEM)

typedef mty::QCD_Model QCD;
REGISTER_MARTY_MODEL("qcd", QCD)

typedef mty::QCD_Model QED;
REGISTER_MARTY_MODEL("qed", QED)

typedef mty::SM_Model StandardModel;
REGISTER_MARTY_MODEL("sm", StandardModel)
