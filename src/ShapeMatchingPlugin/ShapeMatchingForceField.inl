/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
#pragma once

#include <sofa/core/behavior/ForceField.inl>
#include <ShapeMatchingPlugin/ShapeMatchingForceField.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/rmath.h>
#include <assert.h>
#include <iostream>
#include <sofa/simulation/Node.h>

namespace sofa::component::forcefield
{

template<class DataTypes>
ShapeMatchingForceField<DataTypes>::ShapeMatchingForceField()
: l_rotationFinder(initLink("rotationFinder", "link to the rotation finder"))
, d_stiffness(initData(&d_stiffness, (Real)500, "stiffness", "force stiffness"))
{
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);

    Inherit::init();

    if (!l_rotationFinder.get())
    {
        sofa::core::sptr< container::ShapeMatchingRotationFinder<DataTypes>> rotationFinder;
        this->getContext()->get(rotationFinder);
        if (!rotationFinder)
        {
            msg_error() << "did not found any Rotation Finder.";
            return;
        }
        else
        {
            l_rotationFinder.set(rotationFinder);
        }
    }

    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::addForce(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& /* v */)
{
    if (d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
        return;

    sofa::helper::WriteAccessor< core::objectmodel::Data< VecDeriv > > f1 = f;
    sofa::helper::ReadAccessor< core::objectmodel::Data< VecCoord > > p1 = x;

    type::Mat<3,1,Real> q;
    const type::vector<Mat3x3>& vR = l_rotationFinder->getRotations();
    const VecCoord& vcm = l_rotationFinder->getCM();
    const VecCoord& vcm0 = l_rotationFinder->getCM0();
    const typename container::ShapeMatchingRotationFinder<DataTypes>::VecNeighborhood& vNeighborhood = l_rotationFinder->getNeighborhood();
    const VecCoord& x0 = this->mstate->read(sofa::core::ConstVecCoordId::restPosition())->getValue();
    Coord gi, xi;

    f1.resize(p1.size());

    for (unsigned int i=0; i<vNeighborhood.size(); ++i)
    {
        Mat3x3 R = vR[i];
        Coord Xcm = vcm[i];
        Coord Xcm0 = vcm0[i];
        typename container::ShapeMatchingRotationFinder<DataTypes>::Neighborhood::const_iterator it, itEnd;
        for (it = vNeighborhood[i].begin(), itEnd = vNeighborhood[i].end(); it != itEnd ; ++it)
        {
            unsigned int neighbor_i = *it;
            gi = Xcm + R*(x0[neighbor_i] - Xcm0);
            xi = p1[neighbor_i];
            f1[neighbor_i] += (gi-xi)*d_stiffness.getValue();
        }
    }
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx)
{
    if (d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
        return;

    sofa::helper::WriteAccessor< core::objectmodel::Data< VecDeriv > > df1 = df;
    sofa::helper::ReadAccessor< core::objectmodel::Data< VecDeriv > > dp1 = dx;

    const type::vector<Mat3x3>& vDR = l_rotationFinder->getDRotations(dx.getValue());
    const VecCoord& vcm0 = l_rotationFinder->getCM0();
    const typename container::ShapeMatchingRotationFinder<DataTypes>::VecNeighborhood& vNeighborhood = l_rotationFinder->getNeighborhood();
    const VecCoord& x0 = this->mstate->read(sofa::core::ConstVecCoordId::restPosition())->getValue();

    df1.resize(dp1.size());
    const Real fact = d_stiffness.getValue()*mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    for (unsigned int i=0; i<vNeighborhood.size(); ++i)
    {
        const Mat3x3& dR = vDR[i];
        const Coord& Xcm0 = vcm0[i];
        Coord dxcm;

        for(const auto& ni : vNeighborhood[i])
        {
            dxcm += dp1[ni];
        }
        dxcm /= vNeighborhood[i].size();
        for (const auto& ni : vNeighborhood[i])
        {
            Coord dgi = dxcm + dR * (x0[ni] - Xcm0); // + R*(x0[ni] - Xcm0);
            Coord dxi = dp1[ni];
            df1[ni] += (dgi - dxi)*fact;
        }
    }
}

template<class DataTypes>
SReal ShapeMatchingForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams* /*mparams*/, const DataVecCoord&  /*x */) const
{
    // NOT IMPLEMENTED
    return 0.0;
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::draw(const core::visual::VisualParams* )
{
}


} // namespace sofa::component::forcefield
