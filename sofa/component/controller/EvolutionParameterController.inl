/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
//
// C++ Models: EvolutionParameterController
//
// Description:
//
//
// Author: Pierre-Jean Bensoussan, Digital Trainers (2008)
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_CONTROLLER_EVOLUTIONPARAMETERCONTROLLER_INL
#define SOFA_COMPONENT_CONTROLLER_EVOLUTIONPARAMETERCONTROLLER_INL

#include <sofa/component/controller/EvolutionParameterController.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Quat.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/UpdateMappingVisitor.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <SofaBaseTopology/TetrahedronSetGeometryAlgorithms.h>
namespace sofa
{

namespace component
{

namespace controller
{

template <class DataTypes>
EvolutionParameterController<DataTypes>::EvolutionParameterController()
    : f_startDTAppl(initData(&f_startDTAppl, 0, "startDTAppl", "time step when the forcefield starts to be applied"))
    , f_numDTAppl(initData(&f_numDTAppl,0,"numDTAppl", "number of time steps to apply the forcefield"))
    , f_controlFactor(initData(&f_controlFactor,Real(0.0),"factor", "output control factor"))
    , f_compControlFactor(initData(&f_compControlFactor,Real(0.0),"complementFactor", "output control factor substracted from 1.0"))
{
    simulationStep = 0;
}

template <class DataTypes>
void EvolutionParameterController<DataTypes>::init()
{            
}


template <class DataTypes>
void EvolutionParameterController<DataTypes>::applyController(void)
{
}

template <class DataTypes>
void EvolutionParameterController<DataTypes>::onBeginAnimationStep(const double /*dt*/)
{
    if (f_numDTAppl.getValue() == 0 || f_startDTAppl.getValue() == 0) {
        f_controlFactor.setValue(1.0);
        f_compControlFactor.setValue(0.0);
        return;
    }

    Real ff;
    if (simulationStep > f_startDTAppl.getValue()) {
        ff = Real(simulationStep - f_startDTAppl.getValue())/Real(abs(f_numDTAppl.getValue()));
        ff = (ff > 1.0) ? 1.0 : ff;
    }
    else
        ff = 0.0;

    if (f_numDTAppl.getValue() < 0)
        ff = 1.0 - ff;

    f_controlFactor.setValue(ff);
    f_compControlFactor.setValue(1.0-ff);
    simulationStep++;
    sout << this->getName() << "[" << simulationStep << "]  applyFactor = " << f_controlFactor.getValue() << sendl;
}

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_EVOLUTIONPARAMETERCONTROLLER_INL
